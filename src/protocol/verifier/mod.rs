use std::marker::PhantomData;

use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix},
    protocol::{
        Proof, compute_mles, preprocess, sample_random_challenge,
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
        // TODO: separate preprocessing for prover and verifier
        let (m_rq, m_polyring, x_polyring, m_mle, num_vars) = preprocess::<D, F>(m, x);

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

        let x_alpha_eval = mles[2]
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        let ell_eval = mles[3]
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
