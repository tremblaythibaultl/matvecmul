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
use ark_ff::{FftField, Field};
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
        proof: &Proof<D, F>,
    ) -> Result<Vec<F>, ()> {
        let start = Instant::now();
        let mut x_bytes_to_absorb = vec![];
        for ct in x.iter() {
            // consider body only (compressed ciphertext)
            let el = ct.get_ring_element(1).unwrap();
            el.serialize(&mut x_bytes_to_absorb).unwrap();
        }
        let after_process_x = start.elapsed();
        println!("Processing x to absorb took: {:?}", after_process_x);
        println!("\n\n--- VERIFIER TIMINGS ---");

        let start = Instant::now();
        transcript.absorb_bytes_par(&x_bytes_to_absorb);
        let after_absorb = start.elapsed();
        println!("time to absorb x: {:?}", after_absorb);

        let mut beta = transcript.squeeze(1);
        beta.push(beta[0] * beta[0]);

        // assumes RLWE rank is 1
        let start = Instant::now();
        let y_batched = proof
            .y
            .iter()
            .map(|ct| {
                ct.get_ring_element(0)
                    .unwrap()
                    .coeffs
                    .iter()
                    .zip(ct.get_ring_element(1).unwrap().coeffs.iter())
                    .map(|(c0, c1)| {
                        beta[0].mul_by_base_prime_field(&c0) + beta[1].mul_by_base_prime_field(&c1)
                    })
                    .collect::<Vec<F>>()
            })
            .flatten()
            .collect::<Vec<F>>();
        let after_batching_y = start.elapsed();
        println!("time to batch y: {:?}", after_batching_y);

        let start = Instant::now();
        let mut y_bytes_to_absorb = vec![];
        for y_batched_i in y_batched.iter() {
            y_batched_i
                .serialize_uncompressed(&mut y_bytes_to_absorb)
                .unwrap();
        }
        let processed_y = start.elapsed();
        println!("processing batched y took: {:?}", processed_y);

        let start = Instant::now();
        transcript.absorb_bytes_par(&y_bytes_to_absorb);
        let after_absorb_y = start.elapsed();
        println!("time to absorb y : {:?}", after_absorb_y);

        // TODO: split commit/prove in WHIR and absorb commitment to r too

        let m = m_rq.height().ilog2() as usize;
        let start = Instant::now();
        let mut challenges = transcript.squeeze(m + 1);
        let alpha = challenges.pop().unwrap();
        let vec_tau = challenges;
        let after_squeeze = start.elapsed();
        println!("Squeezing challenges took: {:?}", after_squeeze);

        let start = Instant::now();
        // compute unpadded MLEs
        let unpadded_z1_mles = compute_z1_mles(&m_rq, x, z1_num_vars, &beta, &alpha, &vec_tau);
        let after_compute_z1_mles = start.elapsed();
        println!(
            "Computing unpadded Z1 MLEs took: {:?}",
            after_compute_z1_mles
        );

        let powers_of_alpha = unpadded_z1_mles[2].evals();

        let start = Instant::now();
        let (z1_original_claim, z1_final_claim, z1_challenges) = proof
            .z1_sumcheck_proof
            .verify(z1_num_vars, 4, transcript)
            .unwrap();
        let after_verify_z1 = start.elapsed();
        println!("Verifying Z1 took: {:?}", after_verify_z1);

        let start = Instant::now();
        let eq_eval = evaluate_at_challenges(&unpadded_z1_mles[0], &z1_challenges[0..m]);
        let after_eq_eval = start.elapsed();
        println!("Evaluating eq_tau took: {:?}", after_eq_eval);

        let start = Instant::now();
        let x_alpha_eval = evaluate_at_challenges(
            &unpadded_z1_mles[1],
            &z1_challenges[m_rq.height().ilog2() as usize
                ..m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize]
                .to_vec(),
        );
        let after_x_alpha_eval = start.elapsed();
        println!("Evaluating x_alpha took: {:?}", after_x_alpha_eval);

        let start = Instant::now();
        let z1_ell_eval = evaluate_at_challenges(
            &unpadded_z1_mles[2],
            &z1_challenges[m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize..],
        );
        let after_ell_eval = start.elapsed();
        println!("Evaluating ell took: {:?}", after_ell_eval);

        let mut rng = get_rng();

        let whir = Whir::<F>::new(z1_num_vars, &mut rng);
        let start = Instant::now();
        whir.verify(
            &proof.m_mle_proof.as_ref().expect("PCS not implemented"),
            &z1_challenges,
        )
        .expect("m_mle proof does not verify");
        let m_eval = proof.m_mle_proof.as_ref().unwrap().claim;
        let after_m_eval = start.elapsed();
        println!("WHIR verifying m took: {:?}", after_m_eval);

        let z1_eval_at_random_point = eq_eval * m_eval * x_alpha_eval * z1_ell_eval;

        if z1_eval_at_random_point != z1_final_claim {
            return Err(());
        }

        // z_3 verification
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let start = Instant::now();
        let (z3_original_claim, z3_final_claim, z3_challenges) = proof
            .z3_sumcheck_proof
            .verify(z3_num_vars, 3, transcript)
            .unwrap();
        let after_verify_z3 = start.elapsed();
        println!("Verifying Z3 took: {:?}", after_verify_z3);

        let start = Instant::now();
        let z3_ell_eval = evaluate_at_challenges(
            &unpadded_z1_mles[2],
            &z3_challenges[m_rq.height().ilog2() as usize..].to_vec(),
        );
        let after_ell_eval_z3 = start.elapsed();
        println!("Evaluating ell for Z3 took: {:?}", after_ell_eval_z3);

        let r0_eval = proof.r0_mle_proof.as_ref().unwrap().claim;
        let r1_eval = proof.r1_mle_proof.as_ref().unwrap().claim;

        let start = Instant::now();
        let whir = Whir::<F>::new(z3_num_vars, &mut rng);
        let after_whir_setup = start.elapsed();
        println!("Setting up Whir took: {:?}", after_whir_setup);
        let start = Instant::now();
        whir.verify(&proof.r0_mle_proof.as_ref().unwrap(), &z3_challenges)
            .expect("r0_mle proof does not verify");
        whir.verify(&proof.r1_mle_proof.as_ref().unwrap(), &z3_challenges)
            .expect("r1_mle proof does not verify");
        let after_whir_verify = start.elapsed();
        println!("WHIR verifying r0 and r1 took: {:?}", after_whir_verify);

        // eq_tau eval
        let start = Instant::now();
        let eq_eval = evaluate_at_challenges(&unpadded_z1_mles[0], &z3_challenges[0..m]);
        let after_eq_eval_z3 = start.elapsed();
        println!("Evaluating eq_tau for Z3 took: {:?}", after_eq_eval_z3);

        let r_eval = beta[0].mul(&r0_eval) + beta[1].mul(&r1_eval);

        let z3_eval_at_random_point = eq_eval * r_eval * z3_ell_eval;

        if z3_eval_at_random_point != z3_final_claim {
            return Err(());
        }

        // z_2
        let z2_num_vars = m_rq.height().ilog2() as usize;

        let alpha_n_plus_one = (alpha * powers_of_alpha[D - 1]) + F::ONE;
        let alpha_factor = alpha_n_plus_one;

        let start = Instant::now();
        let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << z2_num_vars);
        build_eq_poly(&vec_tau[0..z2_num_vars], &mut eq_tau_mle_evals);
        let after_build_eq = start.elapsed();
        println!("Building eq_tau for Z2 took: {:?}", after_build_eq);

        let start = Instant::now();

        let y_alpha_mle_evals = y_batched
            .chunks(D)
            .map(|chunk| {
                chunk
                    .iter()
                    .zip(powers_of_alpha.iter())
                    .map(|(c, a)| c.mul(a))
                    .sum::<F>()
            })
            .collect::<Vec<F>>();

        let after_y_alpha = start.elapsed();
        println!("Computing y_alpha evaluations took: {:?}", after_y_alpha);

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
    beta: &Vec<F>,
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

    let x_alpha_mle_evals = x
        .par_iter()
        .map(|ct| {
            ct.get_ring_element(0)
                .unwrap()
                .coeffs
                .iter()
                .zip(ct.get_ring_element(1).unwrap().coeffs.iter())
                .map(|(c0, c1)| {
                    beta[0].mul_by_base_prime_field(c0) + beta[1].mul_by_base_prime_field(c1)
                })
                .zip(powers_of_alpha.iter())
                .map(|(c, a)| c.mul(a))
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
// we exploit this structure to evaluate the padded x_alpha MLE at the point defined by `challenges_subvec`
// see MultilinearPolynomial::test_fast_eval_structured for more details
fn evaluate_at_challenges<F: Field>(
    unpadded_mle: &MultilinearPolynomial<F>,
    challenges_subvec: &[F],
) -> F {
    assert_eq!(1 << challenges_subvec.len(), unpadded_mle.evals().len());
    let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << challenges_subvec.len());
    build_eq_poly(
        &challenges_subvec.iter().rev().cloned().collect::<Vec<F>>(),
        &mut eq_x_challenges_mle_evals,
    );

    unpadded_mle
        .evals()
        .iter()
        .zip(eq_x_challenges_mle_evals.iter())
        .map(|(a, b)| *a * *b)
        .sum::<F>()
}
