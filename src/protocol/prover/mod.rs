use std::{marker::PhantomData, time::Instant};

use ark_ff::{FftField, Field};
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    slice::ParallelSlice,
};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof,
        pcs::whir::Whir,
        sumcheck::{multilinear::MultilinearPolynomial, prove},
        transcript::Blake3Transcript,
        utils::{build_eq_poly, sum_over_boolean_hypercube},
    },
    rand::get_rng,
    rlwe::RLWE,
};

pub struct Prover<const D: usize, F: FftField> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F> Prover<D, F>
where
    F: FftField,
{
    pub fn preprocess(
        m: &Matrix<F::BasePrimeField>,
    ) -> (
        Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        Matrix<PolynomialRing<D, F::BasePrimeField>>,
        MultilinearPolynomial<F>,
        MultilinearPolynomial<F::BasePrimeField>,
        usize,
        Blake3Transcript<F>,
    ) {
        // interpret each row as a ring elements and rearrange columns
        let m_rq = m.lift_to_rq::<D>();

        // lift matrix entries to polynomial ring elements of max degree 2*D
        // this performs a clone of the coefficients
        // should look into optimizing
        let m_polyring = m_rq.lift_to_polynomial_ring();

        // this also clones the coefficients. should look into optimizing.
        let (m_mle_evals, num_vars): (Vec<F::BasePrimeField>, usize) = m_rq.to_mle_evals();

        let m_mle_over_base_f = MultilinearPolynomial::new(m_mle_evals.clone(), num_vars);
        let m_mle = MultilinearPolynomial::new(
            m_mle_evals
                .iter()
                .map(|e| F::from_base_prime_field(*e))
                .collect(),
            num_vars,
        );

        let mut transcript = Blake3Transcript::<F>::new();

        for elem in m.data.iter() {
            transcript.absorb(&F::from_base_prime_field(*elem))
        }

        (
            m_rq,
            m_polyring,
            m_mle,
            m_mle_over_base_f,
            num_vars,
            transcript,
        )
    }

    pub fn prove(
        m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        m_polyring: &Matrix<PolynomialRing<D, F::BasePrimeField>>,
        m_mle: &MultilinearPolynomial<F>,
        m_mle_over_base_f: &MultilinearPolynomial<F::BasePrimeField>,
        z1_num_vars: usize,
        transcript: &mut Blake3Transcript<F>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
        include_pcs: bool,
    ) -> Proof<D, F> {
        let mut x_bytes_to_absorb = vec![];
        for ct in x.iter() {
            // consider body only for transcript absorbtion (compressed ciphertext)
            let el = ct.get_ring_element(1).unwrap();
            el.serialize(&mut x_bytes_to_absorb).unwrap();
        }

        transcript.absorb_bytes_par(&x_bytes_to_absorb);

        let mut beta = transcript.squeeze(1);
        beta.push(beta[0] * beta[0]);

        // let start = Instant::now();
        let x_polyring = x
            .iter()
            .map(|elem| elem.lift_to_polynomial_ring())
            .collect::<Vec<_>>();
        // let after_lifting_x = start.elapsed();
        // println!("Lifting x to polynomial ring took: {:?}", after_lifting_x);

        // let start = Instant::now();
        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);
        // let after_mat_vec_mul = start.elapsed();
        // println!(
        //     "Matrix-vector multiplication in polynomial ring took: {:?}",
        //     after_mat_vec_mul
        // );

        // let start = Instant::now();
        let (vec_quotients, _vec_remainders): (
            Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>,
            Vec<Vec<CyclotomicRing<D, F::BasePrimeField>>>,
        ) = y_polyring
            .iter()
            .map(|ct| ct.long_division_by_cyclotomic())
            .collect();

        // let after_division = start.elapsed();
        // println!("time taken for long division: {:?}", after_division);

        let mut rng = get_rng();
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        // let start = Instant::now();
        let (r0_mle_base_prime_f_evals, r1_mle_base_prime_f_evals) =
            compute_r_bpf_mle_evals::<D, F>(&vec_quotients);

        // let after_r_mle_evals = start.elapsed();
        // println!("Computing r mle evals took: {:?}", after_r_mle_evals);

        // let start = Instant::now();
        let r0_mle_over_base_prime_f =
            MultilinearPolynomial::new(r0_mle_base_prime_f_evals.clone(), z3_num_vars);
        let whir_r0_mle = Whir::<F>::new(r0_mle_over_base_prime_f.num_variables(), &mut rng);

        let r1_mle_over_base_prime_f =
            MultilinearPolynomial::new(r1_mle_base_prime_f_evals.clone(), z3_num_vars);
        let whir_r1_mle = Whir::<F>::new(r1_mle_over_base_prime_f.num_variables(), &mut rng);

        let whir_r0_mle_commitment_and_prover_state = if include_pcs {
            Some(whir_r0_mle.commit(&r0_mle_over_base_prime_f))
        } else {
            None
        };

        let whir_r1_mle_commitment_and_prover_state = if include_pcs {
            Some(whir_r1_mle.commit(&r1_mle_over_base_prime_f))
        } else {
            None
        };
        // let after_r_mle_whir_commitment = start.elapsed();
        // println!(
        //     "Committing to r mle with WHIR took: {:?}",
        //     after_r_mle_whir_commitment
        // );
        // let start = Instant::now();
        let whir_m_mle = Whir::<F>::new(m_mle_over_base_f.num_variables(), &mut rng);

        let whir_m_mle_commitment_and_prover_state = if include_pcs {
            Some(whir_m_mle.commit(&m_mle_over_base_f))
        } else {
            None
        };
        // let after_m_mle_whir_commitment = start.elapsed();
        // println!(
        //     "Committing to m mle with WHIR took: {:?}",
        //     after_m_mle_whir_commitment
        // );

        // mat vec mul -- this value should coincide with vec_remainders
        // Maybe it would be better to compute it from `y_polyring`

        // let start = Instant::now();
        let y = m_rq.mat_rlwe_vec_mul(&x);
        // let after_mat_vec_mul_rq = start.elapsed();
        // println!(
        //     "Matrix-vector multiplication in cyclotomic ring took: {:?}",
        //     after_mat_vec_mul_rq
        // );

        // assumes RLWE rank is 1
        // let start = Instant::now();
        let y_batched = y
            .par_iter()
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

        // let after_batching_y = start.elapsed();
        // println!("Batching y took: {:?}", after_batching_y);

        // let start = Instant::now();
        let mut y_bytes_to_absorb = vec![];
        for y_batched_i in y_batched.iter() {
            y_batched_i
                .serialize_uncompressed(&mut y_bytes_to_absorb)
                .unwrap();
        }
        // let after_serializing_y = start.elapsed();
        // println!("Serializing y took: {:?}", after_serializing_y);

        // let start = Instant::now();
        transcript.absorb_bytes_par(&y_bytes_to_absorb);
        // let after_absorb_y = start.elapsed();
        // println!("Absorbing y took: {:?}", after_absorb_y);

        // TODO: absorb the commitment to r too
        // let start = Instant::now();
        let mut challenges = transcript.squeeze(m_rq.height().ilog2() as usize + 1);
        // let after_squeezing_challenges = start.elapsed();
        // println!(
        //     "Squeezing challenges took: {:?}",
        //     after_squeezing_challenges
        // );

        let alpha = challenges.pop().unwrap();
        let vec_tau = challenges;

        // let start = Instant::now();
        let powers_of_alpha: Vec<F> = (0..D)
            .scan(F::one(), |state, _| {
                let result = *state;
                *state *= alpha;
                Some(result)
            })
            .collect();

        // let after_computing_powers_of_alpha = start.elapsed();
        // println!(
        //     "Computing powers of alpha took: {:?}",
        //     after_computing_powers_of_alpha
        // );

        // let start = Instant::now();
        let mut z1_mles = compute_z1_mles(
            &m_rq,
            &m_mle,
            x,
            z1_num_vars,
            &beta,
            &powers_of_alpha,
            &vec_tau,
        );
        // let after_computing_z1_mles = start.elapsed();
        // println!("Computing z1 mles took: {:?}", after_computing_z1_mles);

        // let start = Instant::now();
        let z1_claim = sum_over_boolean_hypercube(&z1_mles);
        // let after_computing_z1_claim = start.elapsed();
        // println!("Computing z1 claim took: {:?}", after_computing_z1_claim);

        // let start = Instant::now();
        let (z1_sumcheck_proof, z1_challenges) =
            prove(z1_claim, &mut z1_mles, z1_num_vars, transcript);
        // let after_proving_z1 = start.elapsed();
        // println!("Proving z1 took: {:?}", after_proving_z1);

        // z_3
        // let start = Instant::now();
        let mut z3_mles = compute_z3_mles::<D, F>(
            &powers_of_alpha,
            &r0_mle_base_prime_f_evals,
            &r1_mle_base_prime_f_evals,
            &beta,
            z3_num_vars,
            &vec_tau,
        );
        // let after_computing_z3_mles = start.elapsed();
        // println!("Computing z3 mles took: {:?}", after_computing_z3_mles);

        // let start = Instant::now();
        let z3_claim = sum_over_boolean_hypercube(&z3_mles);
        // let after_computing_z3_claim = start.elapsed();
        // println!("Computing z3 claim took: {:?}", after_computing_z3_claim);

        // let start = Instant::now();
        let (z3_sumcheck_proof, z3_challenges) =
            prove(z3_claim, &mut z3_mles, z3_num_vars, transcript);
        // let after_proving_z3 = start.elapsed();
        // println!("Proving z3 took: {:?}", after_proving_z3);

        // Whir proofs for r_mle and m_mle.
        let (m_mle_proof, r0_mle_proof, r1_mle_proof) = if include_pcs {
            // let start = Instant::now();
            let (r0_mle_commitment, r0_mle_prover_state) =
                whir_r0_mle_commitment_and_prover_state.unwrap();

            let r0_mle_proof =
                whir_r0_mle.prove(r0_mle_commitment, r0_mle_prover_state, &z3_challenges);
            // let after_r0_mle_proof = start.elapsed();
            // println!(
            //     "time taken to compute r0_mle proof: {:?}",
            //     after_r0_mle_proof
            // );
            // let start = Instant::now();
            let (r1_mle_commitment, r1_mle_prover_state) =
                whir_r1_mle_commitment_and_prover_state.unwrap();

            let r1_mle_proof =
                whir_r1_mle.prove(r1_mle_commitment, r1_mle_prover_state, &z3_challenges);

            // let after_r1_mle_proof = start.elapsed();
            // println!(
            //     "time taken to compute r1_mle proof: {:?}",
            //     after_r1_mle_proof
            // );

            // let start = Instant::now();
            let (m_mle_commitment, m_mle_prover_state) =
                whir_m_mle_commitment_and_prover_state.unwrap();
            let m_mle_proof =
                whir_m_mle.prove(m_mle_commitment, m_mle_prover_state, &z1_challenges);
            // let after_m_mle_proof = start.elapsed();
            // println!("time taken to compute m_mle proof: {:?}", after_m_mle_proof);

            (Some(m_mle_proof), Some(r0_mle_proof), Some(r1_mle_proof))
        } else {
            (None, None, None)
        };

        Proof {
            y,
            z1_sumcheck_proof,
            z3_sumcheck_proof,
            r0_mle_proof,
            r1_mle_proof,
            m_mle_proof,
        }
    }
}

fn compute_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    m_mle: &MultilinearPolynomial<F>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    beta: &Vec<F>,
    powers_of_alpha: &Vec<F>,
    vec_tau: &Vec<F>,
) -> Vec<MultilinearPolynomial<F>> {
    // this might not be very efficient memory-wise, but we need it to keep genericity in the sumcheck prover.
    // TODO: send minimal information to the sumcheck prover and a function describing how to pad the MLEs
    let padded_alpha_mle_evals = (0..m_rq.width() * m_rq.height())
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();

    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, num_vars);

    let x_batched = x
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
                .collect::<Vec<F>>()
        })
        .flatten()
        .collect::<Vec<F>>();

    let x_alpha_mle_evals = x_batched
        .chunks(D)
        .map(|chunk| {
            chunk
                .iter()
                .zip(powers_of_alpha.iter())
                .map(|(c, a)| c.mul(a))
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

fn compute_r_bpf_mle_evals<const D: usize, F: Field>(
    vec_quotients: &Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>,
) -> (Vec<F::BasePrimeField>, Vec<F::BasePrimeField>) {
    let r0_mle_base_prime_f_evals = vec_quotients
        .iter()
        .flat_map(|quotients| &quotients[0].coeffs[0..D]) // only consider the lower order coefficients
        .copied()
        .collect::<Vec<F::BasePrimeField>>();

    let r1_mle_base_prime_f_evals = vec_quotients
        .iter()
        .flat_map(|quotients| &quotients[1].coeffs[0..D]) // only consider the lower order coefficients
        .copied()
        .collect::<Vec<F::BasePrimeField>>();

    (r0_mle_base_prime_f_evals, r1_mle_base_prime_f_evals)
}

fn compute_z3_mles<const D: usize, F: Field>(
    powers_of_alpha: &Vec<F>,
    r0_mle_base_prime_f_evals: &Vec<F::BasePrimeField>,
    r1_mle_base_prime_f_evals: &Vec<F::BasePrimeField>,
    beta: &Vec<F>,
    z3_num_vars: usize,
    vec_tau: &Vec<F>,
) -> Vec<MultilinearPolynomial<F>> {
    let batched_r_mle_evals = r0_mle_base_prime_f_evals
        .iter()
        .zip(r1_mle_base_prime_f_evals.iter())
        .map(|(e0, e1)| beta[0].mul_by_base_prime_field(&e0) + beta[1].mul_by_base_prime_field(&e1))
        .collect::<Vec<F>>();

    let batched_r_mle_evals_len = batched_r_mle_evals.len();
    let batched_r_mle = MultilinearPolynomial::new(batched_r_mle_evals, z3_num_vars);

    let padded_alpha_mle_evals = (0..batched_r_mle_evals_len / D)
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();
    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, z3_num_vars);

    // compute eq polynomial
    let m = (batched_r_mle_evals_len / D).ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; D])
        .flatten()
        .collect::<Vec<F>>();

    let eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, z3_num_vars);

    vec![eq_tau_mle, batched_r_mle, alpha_mle]
}
