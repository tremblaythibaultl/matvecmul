use std::{marker::PhantomData, time::Instant};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix},
    protocol::{
        Proof, pcs::whir::Whir, sumcheck::multilinear::MultilinearPolynomial,
        transcript::Blake3Transcript, utils::build_eq_poly,
    },
    rand::get_rng,
    rlwe::RLWE,
};
use ark_ff::{BigInteger, FftField, Field, PrimeField};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

pub struct Verifier<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F> Verifier<D, F>
where
    F: FftField,
{
    pub fn preprocess(
        m: &Matrix<F::BasePrimeField>,
    ) -> (
        Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        usize,
        Blake3Transcript<F>,
    ) {
        // interpret each row as a ring elements and rearrange columns
        let m_rq = m.lift_to_rq::<D>();

        let num_vars = m_rq.mle_evals_num_vars();

        let mut transcript = Blake3Transcript::<F>::new();

        for elem in m.data.iter() {
            transcript.absorb(&F::from_base_prime_field(*elem))
        }

        (m_rq, num_vars, transcript)
    }

    pub fn verify(
        m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        z1_num_vars: usize,
        transcript: &mut Blake3Transcript<F>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
        proof: Proof<D, F>,
    ) -> Result<Vec<F>, ()> {
        // maybe one could process x's owned data once to obtain a long Vec<u8>
        // and then pass it as reference to the transcript
        // and to the `compute_z1_mles` function

        let start = Instant::now();
        // TODO: look at how to get ring elements without cloning?? not even sure if the clones are optimized out by the compiler.
        let x_bytes_to_absorb = x
            .par_iter()
            .map(|ct| {
                // consider mask only for now
                ct.get_ring_elements()[0]
                    .coeffs
                    .iter()
                    .flat_map(|c| c.into_bigint().to_bytes_le())
                    .collect::<Vec<u8>>()
            })
            .flatten()
            .collect::<Vec<u8>>();
        let after_processing_x = start.elapsed();
        println!("processing x took: {:?}", after_processing_x);

        let start = Instant::now();
        let y_bytes_to_absorb = proof
            .y
            .par_iter()
            .map(|ct| {
                // consider mask only for now
                ct.get_ring_elements()[0]
                    .coeffs
                    .iter()
                    .flat_map(|c| c.into_bigint().to_bytes_le())
                    .collect::<Vec<u8>>()
            })
            .flatten()
            .collect::<Vec<u8>>();
        let after_processing_y = start.elapsed();
        println!("processing y took: {:?}", after_processing_y);

        let start = Instant::now();
        let bytes_to_absorb = &[x_bytes_to_absorb, y_bytes_to_absorb].concat();
        let after_concat = start.elapsed();
        println!("time to concat: {:?}", after_concat);

        let start = Instant::now();
        transcript.absorb_bytes(bytes_to_absorb);
        let after_absorb = start.elapsed();
        println!("time to absorb: {:?}", after_absorb);

        let m = m_rq.height().ilog2() as usize;
        let start = Instant::now();
        let mut challenges = transcript.squeeze(m + 1);
        let alpha = challenges.pop().unwrap();
        let vec_tau = challenges;

        let after_squeeze = start.elapsed();
        println!("Squeezing challenges took: {:?}", after_squeeze);

        let start = Instant::now();
        // compute unpadded MLEs
        let unpadded_z1_mles = compute_z1_mles(&m_rq, x, z1_num_vars, &alpha, &vec_tau);
        let after_compute_z1_mles = start.elapsed();
        println!(
            "Computing unpadded Z1 MLEs took: {:?}",
            after_compute_z1_mles
        );

        let powers_of_alpha = unpadded_z1_mles[2].evals();

        let start = Instant::now();
        let (z1_original_claim, z1_final_claim, z1_challenges) =
            proof.z1_sumcheck_proof.verify(z1_num_vars, 4).unwrap();
        let after_verify_z1 = start.elapsed();
        println!("Verifying Z1 took: {:?}", after_verify_z1);

        let start = Instant::now();
        let eq_eval = evaluate_eq_tau_at_challenges(&unpadded_z1_mles[0], &z1_challenges[0..m]);
        let after_eq_eval = start.elapsed();
        println!("Evaluating eq_tau took: {:?}", after_eq_eval);

        let start = Instant::now();
        let x_alpha_eval = evaluate_x_alpha_at_challenges(
            &unpadded_z1_mles[1],
            &z1_challenges[m_rq.height().ilog2() as usize
                ..m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize]
                .to_vec(),
        );
        let after_x_alpha_eval = start.elapsed();
        println!("Evaluating x_alpha took: {:?}", after_x_alpha_eval);

        let start = Instant::now();
        let z1_ell_eval = evaluate_ell_at_challenges(
            &unpadded_z1_mles[2],
            &z1_challenges[m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize..],
        );
        let after_ell_eval = start.elapsed();
        println!("Evaluating ell took: {:?}", after_ell_eval);

        let mut rng = get_rng();

        let whir = Whir::<F>::new(z1_num_vars, &mut rng);
        let start = Instant::now();
        whir.verify(&proof.m_mle_proof, &z1_challenges)
            .expect("m_mle proof does not verify");
        let m_eval = proof.m_mle_proof.claim;
        let after_m_eval = start.elapsed();
        println!("WHIR verifying m took: {:?}", after_m_eval);

        let z1_eval_at_random_point = eq_eval * m_eval * x_alpha_eval * z1_ell_eval;

        if z1_eval_at_random_point != z1_final_claim {
            return Err(());
        } else {
            println!("z1 Verification successful");
        }

        // z_3 verification
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let start = Instant::now();
        let (z3_original_claim, z3_final_claim, z3_challenges) =
            proof.z3_sumcheck_proof.verify(z3_num_vars, 3).unwrap();
        let after_verify_z3 = start.elapsed();
        println!("Verifying Z3 took: {:?}", after_verify_z3);

        let start = Instant::now();
        let z3_ell_eval = evaluate_ell_at_challenges(
            &unpadded_z1_mles[2],
            &z3_challenges[m_rq.height().ilog2() as usize..].to_vec(),
        );
        let after_ell_eval_z3 = start.elapsed();
        println!("Evaluating ell for Z3 took: {:?}", after_ell_eval_z3);

        let r_eval = proof.r_mle_proof.claim;

        let start = Instant::now();
        let whir = Whir::<F>::new(z3_num_vars, &mut rng);
        let after_whir_setup = start.elapsed();
        println!("Setting up Whir took: {:?}", after_whir_setup);
        let start = Instant::now();
        whir.verify(&proof.r_mle_proof, &z3_challenges)
            .expect("r_mle proof does not verify");
        let after_whir_verify = start.elapsed();
        println!("WHIR Verifying r took: {:?}", after_whir_verify);
        assert_eq!(r_eval, proof.r_mle_proof.claim);

        // eq_tau eval
        let start = Instant::now();
        let eq_eval = evaluate_eq_tau_at_challenges(&unpadded_z1_mles[0], &z1_challenges[0..m]);
        let after_eq_eval_z3 = start.elapsed();
        println!("Evaluating eq_tau for Z3 took: {:?}", after_eq_eval_z3);

        let z3_eval_at_random_point = eq_eval * r_eval * z3_ell_eval;

        if z3_eval_at_random_point != z3_final_claim {
            return Err(());
        } else {
            println!("z3 Verification successful");
        }

        // z_2
        let z2_num_vars = m_rq.height().ilog2() as usize;

        let alpha_n_plus_one = (alpha * powers_of_alpha[D - 1]) + F::ONE;
        let alpha_factor = alpha_n_plus_one; //.pow([1 << z2_num_vars as u64]);

        let start = Instant::now();
        let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << z2_num_vars);
        build_eq_poly(&vec_tau[0..z2_num_vars], &mut eq_tau_mle_evals);
        let after_build_eq = start.elapsed();
        println!("Building eq_tau for Z2 took: {:?}", after_build_eq);

        // focus on the first component of the ciphertext only for now.
        let start = Instant::now();
        let y_alpha_mle_evals = proof
            .y
            .par_iter()
            .map(|ct| {
                ct.get_ring_elements()[1]
                    .coeffs
                    .iter()
                    .zip(powers_of_alpha.iter())
                    .map(|(c, a)| a.mul_by_base_prime_field(c))
                    .sum::<F>()
            })
            .collect::<Vec<F>>();
        let after_y_alpha = start.elapsed();
        println!("Computing y_alpha evaluations took: {:?}", after_y_alpha);

        assert_eq!(y_alpha_mle_evals.len(), eq_tau_mle_evals.len());

        let start = Instant::now();
        let z_2 = y_alpha_mle_evals
            .iter()
            .zip(eq_tau_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();
        let after_z2 = start.elapsed();
        println!("Evaluating Z2 took: {:?}", after_z2);

        assert_eq!(
            z1_original_claim - z_2 - (alpha_factor * z3_original_claim),
            F::ZERO
        );

        Ok(vec![
            z1_original_claim,
            z_2,
            alpha_factor * z3_original_claim,
        ])
    }
}

pub fn compute_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    _num_vars: usize,
    alpha: &F,
    vec_tau: &Vec<F>,
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
            ct.get_ring_elements()[1]
                .coeffs
                .iter()
                .zip(powers_of_alpha.iter())
                .map(|(c, a)| a.mul_by_base_prime_field(c))
                .sum::<F>()
        })
        .collect::<Vec<F>>();

    let x_alpha_mle_num_vars = x_alpha_mle_evals.len().ilog2() as usize;
    let unpadded_x_alpha_mle = MultilinearPolynomial::new(x_alpha_mle_evals, x_alpha_mle_num_vars);

    let m = m_rq.height().ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);
    let unpadded_eq_tau_mle = MultilinearPolynomial::new(eq_tau_mle_evals, m);

    vec![
        unpadded_eq_tau_mle,
        unpadded_x_alpha_mle,
        unpadded_alpha_mle,
    ]
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

fn evaluate_eq_tau_at_challenges<F: Field>(
    unpadded_eq_tau_mle: &MultilinearPolynomial<F>,
    challenges_subvec: &[F],
) -> F {
    let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << challenges_subvec.len());
    build_eq_poly(&challenges_subvec, &mut eq_x_challenges_mle_evals);

    unpadded_eq_tau_mle
        .evals()
        .iter()
        .zip(eq_x_challenges_mle_evals.iter())
        .map(|(a, b)| *a * *b)
        .sum::<F>()
}
