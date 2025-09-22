use std::marker::PhantomData;

use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof, sample_random_challenge,
        sumcheck::{multilinear::MultilinearPolynomial, prove},
        utils::{build_eq_poly, sum_over_boolean_hypercube},
    },
    rlwe::RLWE,
};

pub struct Prover<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: Field> Prover<D, F> {
    pub fn prove(
        m: &Matrix<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    ) -> Proof<D, F> {
        let (m_rq, m_polyring, x_polyring, m_mle, num_vars) = preprocess::<D, F>(m, x);

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        //reduce each polynomial ring element to cyclotomic ring using long division
        let mut vec_quotients =
            Vec::<Vec<PolynomialRing<D, F::BasePrimeField>>>::with_capacity(y_polyring.len());
        let mut vec_remainders =
            Vec::<Vec<CyclotomicRing<D, F::BasePrimeField>>>::with_capacity(y_polyring.len());

        for elem in y_polyring {
            let (quotients, remainders) = elem.long_division_by_cyclotomic();
            vec_quotients.push(quotients);
            vec_remainders.push(remainders);
        }

        // mat vec mul -- this value should coincide with vec_remainders
        let y = m_rq.mat_rlwe_vec_mul(&x);

        println!("y: {:#?}", y);
        println!("vec_remainders: {:#?}", vec_remainders);
        println!("vec_quotients: {:#?}", vec_quotients);

        // TODO: implement fiat shamir
        let alpha = sample_random_challenge::<F>(true);
        let tau = sample_random_challenge::<F>(true);

        let mut mles = compute_mles(&m_rq, &m_mle, x, num_vars, &alpha, &tau);

        let claim = sum_over_boolean_hypercube(&mles);

        let (z1_sumcheck_proof, challenges) = prove(claim, &mut mles, num_vars);

        // We probably will be able to batch the two sumchecks (z_1 and z_3). Not clear how yet.

        // z_3
        // let r_alpha_mle_evals = vec_quotients
        //     .iter()
        //     .map(|elem| {
        //         elem.iter()
        //             .map(|poly| {
        //                 poly.coeffs
        //                     .iter()
        //                     .zip(powers_of_alpha.iter())
        //                     .map(|(c, a)| a.mul_by_base_prime_field(c))
        //                     .sum::<F>()
        //             })
        //             .collect::<Vec<F>>()
        //     })
        //     .flatten()
        //     .collect::<Vec<F>>();

        // let r_alpha_mle: MultilinearPolynomial<F> = MultilinearPolynomial::new(
        //     r_alpha_mle_evals,
        //     (vec_quotients.len() * vec_quotients[0].len())
        //         .next_power_of_two()
        //         .ilog2() as usize,
        // );

        Proof {
            y,
            r: vec_quotients,
            z1_sumcheck_proof,
        }
    }
}

pub fn preprocess<const D: usize, F: Field>(
    m: &Matrix<F::BasePrimeField>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
) -> (
    Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    Matrix<PolynomialRing<D, F::BasePrimeField>>,
    Vec<RLWE<PolynomialRing<D, F::BasePrimeField>>>,
    MultilinearPolynomial<F>,
    usize,
) {
    // interpret each row as a ring elements and rearrange columns
    let m_rq = m.lift_to_rq::<D>();

    // lift matrix entries to polynomial ring elements of max degree 2*D
    // this performs a clone of the coefficients
    // should look into optimizing
    let m_polyring = m_rq.lift_to_polynomial_ring();

    let x_polyring = x
        .iter()
        .map(|elem| elem.lift_to_polynomial_ring())
        .collect::<Vec<_>>();

    // this also clones the coefficients. should look into optimizing.
    let (m_mle_evals, num_vars): (Vec<F::BasePrimeField>, usize) = m_rq.to_mle_evals();

    let m_mle = MultilinearPolynomial::new(
        m_mle_evals
            .into_iter()
            .map(|e| F::from_base_prime_field(e))
            .collect(),
        num_vars,
    );

    (m_rq, m_polyring, x_polyring, m_mle, num_vars)
}

pub fn compute_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    m_mle: &MultilinearPolynomial<F>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    alpha: &F,
    tau: &F,
) -> Vec<MultilinearPolynomial<F>> {
    let powers_of_alpha: Vec<F> = (0..D)
        .scan(F::one(), |state, _| {
            let result = *state;
            *state *= alpha;
            Some(result)
        })
        .collect();

    // this might not be very efficient memory-wise, but we need it to keep genericity in the sumcheck prover.
    // TODO: send minimal information to the sumcheck prover and a function describing how to pad the MLEs
    let padded_alpha_mle_evals = (0..m_rq.width() * m_rq.height())
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();

    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, num_vars);

    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank
    let x_alpha_mle_evals = x
        .iter()
        .map(|ct| {
            ct.get_ring_elements()[0]
                .coeffs
                .iter()
                .zip(powers_of_alpha.iter())
                .map(|(c, a)| a.mul_by_base_prime_field(c))
                .sum::<F>()
        })
        .collect::<Vec<F>>();

    // pad x_alpha_mle_evals
    let padded_x_alpha_mle_evals = (0..m_rq.height() as usize)
        .flat_map(|_| {
            x_alpha_mle_evals
                .iter()
                .map(|e| vec![*e; D])
                .flatten()
                .collect::<Vec<F>>()
        })
        .collect::<Vec<F>>();

    let x_alpha_mle: MultilinearPolynomial<F> =
        MultilinearPolynomial::new(padded_x_alpha_mle_evals, num_vars);

    // compute eq polynomial
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(num_vars);
    let vec_tau = vec![*tau; num_vars];
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let eq_tau_mle = MultilinearPolynomial::new(eq_tau_mle_evals, num_vars);

    vec![
        eq_tau_mle,
        m_mle.clone(),
        x_alpha_mle.clone(),
        alpha_mle.clone(),
    ]
}
