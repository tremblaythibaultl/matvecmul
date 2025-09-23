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

        let unpadded_mles = compute_z1_mles(&m_rq, x, num_vars, &alpha);

        let (original_claim, final_claim, challenges) =
            proof.z1_sumcheck_proof.verify(num_vars, 4).unwrap();

        let eq_eval = eq_eval(&vec_tau, &challenges);

        let x_alpha_eval = evaluate_x_alpha_at_challenges(
            &unpadded_mles[0],
            &challenges[m_rq.height().ilog2() as usize
                ..m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize]
                .to_vec(),
        );

        let ell_eval = evaluate_ell_at_challenges(
            &unpadded_mles[1],
            &challenges[m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize..],
        );

        // This does NOT need to be computed by the verifier
        // The value should come from the PCS.
        // Included for test purposes
        let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << num_vars);
        build_eq_poly(&challenges, &mut eq_x_challenges_mle_evals);

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

pub fn compute_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    alpha: &F,
) -> Vec<MultilinearPolynomial<F>> {
    let powers_of_alpha: Vec<F> = (0..D)
        .scan(F::one(), |state, _| {
            let result = *state;
            *state *= alpha;
            Some(result)
        })
        .collect();

    let powers_of_alpha_num_vars = powers_of_alpha.len().ilog2() as usize;
    let unpadded_alpha_mle =
        MultilinearPolynomial::new(powers_of_alpha.clone(), powers_of_alpha_num_vars);

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

    let x_alpha_mle_num_vars = x_alpha_mle_evals.len().ilog2() as usize;
    let unpadded_x_alpha_mle = MultilinearPolynomial::new(x_alpha_mle_evals, x_alpha_mle_num_vars);

    vec![unpadded_x_alpha_mle, unpadded_alpha_mle]
}

// The padded x_alpha MLE has `num_vars` variables.
// Its evaluations `evals` are the `D` repetitions of each evaluation of the unpadded x_alpha MLE, concatenated `m` times.
// In other words, `evals[0] = evals[1] = ... = evals[D-1], evals[D] = evals[D+1] = ... = evals[2D-1]` up to `evals [t*D-1]`, then
// `evals[0] = evals[t*D]` and so on and so forth
fn evaluate_x_alpha_at_challenges<F: Field>(
    unpadded_x_alpha_mle: &MultilinearPolynomial<F>,
    challenges_subvec: &[F],
) -> F {
    let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << challenges_subvec.len());
    build_eq_poly(&challenges_subvec, &mut eq_x_challenges_mle_evals);

    unpadded_x_alpha_mle
        .evals()
        .iter()
        .zip(eq_x_challenges_mle_evals.iter())
        .map(|(a, b)| *a * *b)
        .sum::<F>()
}

// Similar reasoning as `evaluate_x_alpha_at_challenges`
fn evaluate_ell_at_challenges<F: Field>(
    unpadded_ell_mle: &MultilinearPolynomial<F>,
    challenges_subvec: &[F],
) -> F {
    let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << challenges_subvec.len());
    build_eq_poly(&challenges_subvec, &mut eq_x_challenges_mle_evals);

    unpadded_ell_mle
        .evals()
        .iter()
        .zip(eq_x_challenges_mle_evals.iter())
        .map(|(a, b)| *a * *b)
        .sum::<F>()
}
