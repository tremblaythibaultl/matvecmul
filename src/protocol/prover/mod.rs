use std::marker::PhantomData;

use ark_ff::{FftField, Field};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof,
        prover::whir::Whir,
        sample_random_challenge,
        sumcheck::{multilinear::MultilinearPolynomial, prove},
        utils::{build_eq_poly, sum_over_boolean_hypercube},
    },
    rlwe::RLWE,
};

pub mod whir;

pub struct Prover<const D: usize, F: FftField> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: FftField> Prover<D, F> {
    pub fn prove(
        m: &Matrix<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    ) -> Proof<D, F> {
        let (m_rq, m_polyring, x_polyring, m_mle, z1_num_vars) = preprocess::<D, F>(m, x);

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

        let powers_of_alpha: Vec<F> = (0..D)
            .scan(F::one(), |state, _| {
                let result = *state;
                *state *= alpha;
                Some(result)
            })
            .collect();

        let mut z1_mles = compute_z1_mles(&m_rq, &m_mle, x, z1_num_vars, &powers_of_alpha, &tau);

        let z1_claim = sum_over_boolean_hypercube(&z1_mles);

        let (z1_sumcheck_proof, z1_challenges) = prove(z1_claim, &mut z1_mles, z1_num_vars);

        // We probably will be able to batch the two sumchecks (z_1 and z_3). Not clear how yet.

        // z_3
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let mut z3_mles = compute_z3_mles(&vec_quotients, &powers_of_alpha, z3_num_vars, &tau);

        let z3_claim = sum_over_boolean_hypercube(&z3_mles);

        // included for test purposes only. The real protocol should commit to this using a PCS
        let r_mle = z3_mles[1].clone();

        let (z3_sumcheck_proof, z3_challenges) = prove(z3_claim, &mut z3_mles, z3_num_vars);

        // Whir proof for r_mle.
        let mut rng = ark_std::test_rng();
        let whir = Whir::<F>::new(r_mle.num_variables(), &mut rng);
        let r_mle_proof = whir.prove(&r_mle, &z3_challenges);

        // sanity check with z_2

        // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank
        let y_alpha_mle_evals = y
            .iter()
            .map(|ct| {
                ct.get_ring_elements()[1]
                    .coeffs
                    .iter()
                    .zip(powers_of_alpha.iter())
                    .map(|(c, a)| a.mul_by_base_prime_field(c))
                    .sum::<F>()
            })
            .collect::<Vec<F>>();

        let z2_num_vars = m_rq.height().ilog2() as usize;
        let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << z2_num_vars);
        let vec_tau = vec![tau; z2_num_vars];
        build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

        let z_2 = sum_over_boolean_hypercube(&[
            MultilinearPolynomial::new(y_alpha_mle_evals, z2_num_vars),
            MultilinearPolynomial::new(eq_tau_mle_evals, z2_num_vars),
        ]);

        Proof {
            y,
            r: vec_quotients,
            z1_sumcheck_proof,
            z3_sumcheck_proof,
            r_mle,
            r_mle_proof,
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

fn compute_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    m_mle: &MultilinearPolynomial<F>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    powers_of_alpha: &Vec<F>,
    tau: &F,
) -> Vec<MultilinearPolynomial<F>> {
    // this might not be very efficient memory-wise, but we need it to keep genericity in the sumcheck prover.
    // TODO: send minimal information to the sumcheck prover and a function describing how to pad the MLEs
    let padded_alpha_mle_evals = (0..m_rq.width() * m_rq.height())
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();

    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, num_vars);

    // only consider the body for now. will probably need to run the protocol K+1 times where K is the RLWE rank
    let x_alpha_mle_evals = x
        .iter()
        .map(|ct| {
            ct.get_ring_elements()[1]
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
    let m = m_rq.height().ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    let vec_tau = vec![*tau; m];
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; m_rq.width() as usize * D])
        .flatten()
        .collect::<Vec<F>>();

    let eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, num_vars);

    vec![
        eq_tau_mle,
        m_mle.clone(),
        x_alpha_mle.clone(),
        alpha_mle.clone(),
    ]
}

fn compute_z3_mles<const D: usize, F: Field>(
    vec_quotients: &Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>,
    powers_of_alpha: &Vec<F>,
    z3_num_vars: usize,
    tau: &F,
) -> Vec<MultilinearPolynomial<F>> {
    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank

    let r_mle_evals = vec_quotients
        .iter()
        .map(|quotients| {
            quotients[1].coeffs[0..D] // only consider the lower order coefficients
                .iter()
                .map(|c| F::from_base_prime_field(*c))
                .collect::<Vec<F>>()
        })
        .flatten()
        .collect::<Vec<F>>();
    let r_mle_evals_len = r_mle_evals.len();

    let r_mle = MultilinearPolynomial::new(r_mle_evals, z3_num_vars);

    let padded_alpha_mle_evals = (0..r_mle_evals_len / D)
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();
    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, z3_num_vars);

    // compute eq polynomial
    let m = (r_mle_evals_len / D).ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    let vec_tau = vec![*tau; m];
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; D])
        .flatten()
        .collect::<Vec<F>>();

    let eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, z3_num_vars);

    vec![eq_tau_mle, r_mle, alpha_mle]
}
