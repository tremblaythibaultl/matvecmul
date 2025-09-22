use std::marker::PhantomData;

use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof, sample_random_challenge,
        sumcheck::multilinear::MultilinearPolynomial,
        utils::{build_eq_poly, eq_eval},
    },
    rlwe::RLWE,
};

pub struct Verifier<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: Field> Verifier<D, F> {
    pub fn verify(
        m: &Matrix<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
        proof: Proof<D, F>,
    ) -> Result<Vec<F::BasePrimeField>, ()> {
        let (m_rq, m_mle, num_vars) = preprocess::<D, F>(m, x);

        // TODO: implement fiat shamir
        let alpha = sample_random_challenge::<F>(true);
        let tau = sample_random_challenge::<F>(true);
        let vec_tau = vec![tau; num_vars];

        // compute padded MLEs

        let mles = compute_mles(&m_rq, &m_mle, x, num_vars, &alpha, &tau);

        let (original_claim, final_claim, challenges) =
            proof.z1_sumcheck_proof.verify(num_vars, 4).unwrap();

        let eq_eval = eq_eval(&vec_tau, &challenges);

        let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(num_vars);
        build_eq_poly(&challenges, &mut eq_x_challenges_mle_evals);

        let x_alpha_eval = mles[0]
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        let ell_eval = mles[1]
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        // This does NOT need to be computed by the verifier
        // The value should come from the PCS.
        // Included for test purposes
        let m_eval = m_mle
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        let eval_at_random_point = eq_eval * m_eval * x_alpha_eval * ell_eval;

        if eval_at_random_point != final_claim {
            return Err(());
        } else {
            println!("Verification successful");
        }
        Ok(original_claim.to_base_prime_field_elements().collect())
    }
}

pub fn preprocess<const D: usize, F: Field>(
    m: &Matrix<F::BasePrimeField>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
) -> (
    Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    MultilinearPolynomial<F>,
    usize,
) {
    // interpret each row as a ring elements and rearrange columns
    let m_rq = m.lift_to_rq::<D>();

    // this clones the coefficients. should look into optimizing.
    let (m_mle_evals, num_vars): (Vec<F::BasePrimeField>, usize) = m_rq.to_mle_evals();

    let m_mle = MultilinearPolynomial::new(
        m_mle_evals
            .into_iter()
            .map(|e| F::from_base_prime_field(e))
            .collect(),
        num_vars,
    );

    (m_rq, m_mle, num_vars)
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

    vec![x_alpha_mle.clone(), alpha_mle.clone()]
}
