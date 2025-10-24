use std::{marker::PhantomData, ops::Mul, time::Instant};

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
        proof: Proof<D, F>,
    ) -> Result<Vec<F>, ()> {
        let start = Instant::now();
        let mut x_bytes_to_absorb = vec![];
        for ct in x.iter() {
            // consider mask only for now
            let el = ct.get_ring_element(0).unwrap();
            el.serialize(&mut x_bytes_to_absorb).unwrap();
        }

        let after_processing_x = start.elapsed();
        println!("processing x took: {:?}", after_processing_x);

        let start = Instant::now();
        let mut y_bytes_to_absorb = vec![];
        for ct in proof.y.iter() {
            // consider mask only for now
            let el = ct.get_ring_element(0).unwrap();
            el.serialize(&mut y_bytes_to_absorb).unwrap();
        }
        let after_processing_y = start.elapsed();
        println!("processing y took: {:?}", after_processing_y);

        let start = Instant::now();
        let bytes_to_absorb = [x_bytes_to_absorb, y_bytes_to_absorb].concat();
        let after_concat = start.elapsed();
        println!("time to concat: {:?}", after_concat);

        let start = Instant::now();
        transcript.absorb_bytes_par(&bytes_to_absorb);
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

        let start = Instant::now();
        let mut padded_z1_mles = compute_padded_z1_mles(&m_rq, x, z1_num_vars, &alpha, &vec_tau);
        let after_compute_padded_z1_mles = start.elapsed();
        println!(
            "Computing padded Z1 MLEs took: {:?}",
            after_compute_padded_z1_mles
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
        for chal in z1_challenges.iter() {
            // TODO: par_iter_mut()?
            for mle in padded_z1_mles.iter_mut() {
                mle.bind_to_challenge(&chal);
            }
        }
        let after_bind_challenges = start.elapsed();
        println!("Binding Z1 challenges took: {:?}", after_bind_challenges);

        let eq_eval = padded_z1_mles[0].evals()[0];
        let x_alpha_eval = padded_z1_mles[1].evals()[0];
        let z1_ell_eval = padded_z1_mles[2].evals()[0];

        let mut rng = get_rng();

        let whir = Whir::<F>::new(z1_num_vars, &mut rng);
        let start = Instant::now();
        whir.verify(&proof.m_mle_proof, &z1_challenges)
            .expect("m_mle proof does not verify");
        let m_eval = proof.m_mle_proof.claim;
        let after_m_eval = start.elapsed();
        println!("WHIR verifying m took: {:?}", after_m_eval);

        let z1_eval_at_random_point = eq_eval * m_eval * x_alpha_eval * z1_ell_eval;

        println!(
            "z1_final_claim: {:?}, z1_eval_at_random_point: {:?}",
            z1_final_claim, z1_eval_at_random_point
        );

        if z1_eval_at_random_point != z1_final_claim {
            return Err(());
        } else {
            println!("z1 Verification successful");
        }

        // z_3 verification
        let m = m_rq.height().ilog2() as usize;
        let t = m_rq.width().ilog2() as usize;

        let z3_num_vars = z1_num_vars - t;

        let start = Instant::now();
        let (z3_original_claim, z3_final_claim, z3_challenges) = proof
            .z3_sumcheck_proof
            .verify(z3_num_vars, 3, transcript)
            .unwrap();
        let after_verify_z3 = start.elapsed();
        println!("Verifying Z3 took: {:?}", after_verify_z3);

        let mut padded_z3_mles =
            compute_padded_z3_mles::<D, F>(&powers_of_alpha.to_vec(), &vec_tau, m);

        let start = Instant::now();
        for chal in z3_challenges.iter() {
            // TODO: par_iter_mut()?
            for mle in padded_z3_mles.iter_mut() {
                mle.bind_to_challenge(&chal);
            }
        }
        let after_bind_challenges = start.elapsed();
        println!("Binding Z3 challenges took: {:?}", after_bind_challenges);

        println!("padded_z3_mles[0] len: {}", padded_z3_mles[0].evals().len());
        println!("padded_z3_mles[1] len: {}", padded_z3_mles[1].evals().len());

        let eq_eval = padded_z3_mles[0].evals()[0];
        let z3_ell_eval = padded_z3_mles[1].evals()[0];

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
                ct.get_ring_element(1)
                    .unwrap()
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

pub fn compute_padded_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
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
    // pad alpha_mle_evals
    let padded_alpha_mle_evals = (0..m_rq.width() * m_rq.height())
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();

    let padded_alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, num_vars);

    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank
    let x_alpha_mle_evals = x
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

    let padded_x_alpha_mle_evals = (0..m_rq.height() as usize)
        .flat_map(|_| {
            x_alpha_mle_evals
                .iter()
                .map(|e| vec![*e; D])
                .flatten()
                .collect::<Vec<F>>()
        })
        .collect::<Vec<F>>();

    let padded_x_alpha_mle: MultilinearPolynomial<F> =
        MultilinearPolynomial::new(padded_x_alpha_mle_evals, num_vars);

    let m = m_rq.height().ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);
    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; m_rq.width() as usize * D])
        .flatten()
        .collect::<Vec<F>>();

    let padded_eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, num_vars);
    vec![padded_eq_tau_mle, padded_x_alpha_mle, padded_alpha_mle]
}

fn compute_padded_z3_mles<const D: usize, F: Field>(
    powers_of_alpha: &Vec<F>,
    vec_tau: &Vec<F>,
    m: usize,
) -> Vec<MultilinearPolynomial<F>> {
    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank

    let z3_num_vars = m + D.ilog2() as usize;

    let padded_alpha_mle_evals = (0..1 << m)
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();
    let padded_alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, z3_num_vars);

    // compute eq polynomial
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; D])
        .flatten()
        .collect::<Vec<F>>();

    let padded_eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, z3_num_vars);

    vec![padded_alpha_mle, padded_eq_tau_mle]
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
        .par_iter()
        .map(|ct| {
            ct.get_ring_element(1)
                .unwrap()
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
