use std::{marker::PhantomData, time::Instant};

use ark_ff::{FftField, Field};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof,
        pcs::whir::{Whir, WhirProof},
        sumcheck::{multilinear::MultilinearPolynomial, prove},
        transcript::Blake3Transcript,
        utils::{build_eq_poly, sum_over_boolean_hypercube},
    },
    rand::get_rng,
    rlwe::RLWE,
};

pub struct Z3Mles<F: Field> {
    pub mles_over_f: Vec<MultilinearPolynomial<F>>,
    pub r_mle_over_base_prime_f: MultilinearPolynomial<F::BasePrimeField>,
}

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
        // let start = Instant::now();
        let x_polyring = x
            .iter()
            .map(|elem| elem.lift_to_polynomial_ring())
            .collect::<Vec<_>>();
        // let after_lifting_to_polyring = start.elapsed();
        // println!(
        //     "time taken to lift to polyring: {:?}",
        //     after_lifting_to_polyring
        // );

        // let start = Instant::now();

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        // let after_mat_vec_rlwe_mul = start.elapsed();
        // println!(
        //     "after_mat_vec_rlwe_mul POLYRING: {:?}",
        //     after_mat_vec_rlwe_mul
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
        // println!("after_div: {:?}", after_division);

        let mut rng = get_rng();
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;
        let r_mle_base_prime_f_evals = vec_quotients
            .iter()
            .flat_map(|quotients| &quotients[1].coeffs[0..D]) // only consider the lower order coefficients
            .copied()
            .collect::<Vec<F::BasePrimeField>>();

        let r_mle_over_base_prime_f =
            MultilinearPolynomial::new(r_mle_base_prime_f_evals, z3_num_vars);

        let whir_r_mle = Whir::<F>::new(r_mle_over_base_prime_f.num_variables(), &mut rng);

        let start = Instant::now();

        let whir_r_mle_commitment_and_prover_state = if include_pcs {
            Some(whir_r_mle.commit(&r_mle_over_base_prime_f))
        } else {
            None
        };

        let after_r_commit = start.elapsed();
        println!("time taken to commit r_mle: {:?}", after_r_commit);

        let whir_m_mle = Whir::<F>::new(m_mle_over_base_f.num_variables(), &mut rng);

        let start = Instant::now();
        let whir_m_mle_commitment_and_prover_state = if include_pcs {
            Some(whir_m_mle.commit(&m_mle_over_base_f))
        } else {
            None
        };
        let after_m_commit = start.elapsed();
        println!("time taken to commit m_mle: {:?}", after_m_commit);

        // let start = Instant::now();

        // mat vec mul -- this value should coincide with vec_remainders
        let y = m_rq.mat_rlwe_vec_mul(&x);

        // let after_matvec_mul = start.elapsed();
        // println!("after_mat_vec_mul RLWE: {:?}", after_matvec_mul);

        // let start = Instant::now();
        let mut x_bytes_to_absorb = vec![];
        for ct in x.iter() {
            // consider mask only for now
            let el = ct.get_ring_element(0).unwrap();
            el.serialize(&mut x_bytes_to_absorb).unwrap();
        }
        // let after_processing_x = start.elapsed();
        // println!("processing x took: {:?}", after_processing_x);

        // let start = Instant::now();
        let mut y_bytes_to_absorb = vec![];
        for ct in y.iter() {
            // consider mask only for now
            let el = ct.get_ring_element(0).unwrap();
            el.serialize(&mut y_bytes_to_absorb).unwrap();
        }
        // let after_processing_y = start.elapsed();
        // println!("processing y took: {:?}", after_processing_y);

        // let start = Instant::now();
        let bytes_to_absorb = &[x_bytes_to_absorb, y_bytes_to_absorb].concat();
        // let after_concat = start.elapsed();
        // println!("time to concat: {:?}", after_concat);

        // let start = Instant::now();
        transcript.absorb_bytes_par(&bytes_to_absorb);
        // let after_absorb = start.elapsed();
        // println!("time to absorb: {:?}", after_absorb);

        // TODO: absorb the commitment to r too
        // let start = Instant::now();
        let mut challenges = transcript.squeeze(m_rq.height().ilog2() as usize + 1);
        // let after_squeeze = start.elapsed();
        // println!("Squeezing challenges took: {:?}", after_squeeze);

        let alpha = challenges.pop().unwrap();
        let vec_tau = challenges;

        let powers_of_alpha: Vec<F> = (0..D)
            .scan(F::one(), |state, _| {
                let result = *state;
                *state *= alpha;
                Some(result)
            })
            .collect();

        // let start = Instant::now();
        let mut z1_mles =
            compute_z1_mles(&m_rq, &m_mle, x, z1_num_vars, &powers_of_alpha, &vec_tau);
        // let after_z1_mles = start.elapsed();
        // println!("time to compute z1 mles: {:?}", after_z1_mles);

        // let start = Instant::now();
        let z1_claim = sum_over_boolean_hypercube(&z1_mles);
        // let after_z1_claim = start.elapsed();
        // println!(
        //     "time to compute z1 claim over hypercube: {:?}",
        //     after_z1_claim
        // );

        // let start = Instant::now();
        let (z1_sumcheck_proof, z1_challenges) =
            prove(z1_claim, &mut z1_mles, z1_num_vars, transcript);
        // let after_z1_sc = start.elapsed();
        // println!("time taken for 1st sumcheck: {:?}", after_z1_sc);

        // We probably will be able to batch the two sumchecks (z_1 and z_3). Not clear how yet.

        // z_3

        // let start = Instant::now();
        let mut z3_mles = compute_z3_mles(&vec_quotients, &powers_of_alpha, z3_num_vars, &vec_tau);
        // let after_z3_mles = start.elapsed();
        // println!("time taken to compute z3 mles: {:?}", after_z3_mles);

        // let start = Instant::now();
        let z3_claim = sum_over_boolean_hypercube(&z3_mles.mles_over_f);
        // let after_z3_hypercupe = start.elapsed();
        // println!(
        //     "time taken to compute z3 over hypercube: {:?}",
        //     after_z3_hypercupe
        // );

        // let start = Instant::now();
        let (z3_sumcheck_proof, z3_challenges) =
            prove(z3_claim, &mut z3_mles.mles_over_f, z3_num_vars, transcript);
        // let after_z3_sc = start.elapsed();
        // println!("time taken for 2nd sumcheck: {:?}", after_z3_sc);

        // Whir proofs for r_mle and m_mle.

        let (m_mle_proof, r_mle_proof) = if include_pcs {
            // let start = Instant::now();
            let (r_mle_commitment, r_mle_prover_state) =
                whir_r_mle_commitment_and_prover_state.unwrap();

            let r_mle_proof =
                whir_r_mle.prove(r_mle_commitment, r_mle_prover_state, &z3_challenges);

            // let r_mle_proof = whir_r_mle.commit_and_prove(&z3_mles.r_mle_over_base_prime_f, &z3_challenges);
            // let after_r_mle_proof = start.elapsed();
            // println!("time taken to compute r_mle proof: {:?}", after_r_mle_proof);

            // let start = Instant::now();
            let (m_mle_commitment, m_mle_prover_state) =
                whir_m_mle_commitment_and_prover_state.unwrap();
            let m_mle_proof =
                whir_m_mle.prove(m_mle_commitment, m_mle_prover_state, &z1_challenges);
            // let after_m_mle_proof = start.elapsed();
            // println!("time taken to compute m_mle proof: {:?}", after_m_mle_proof);

            (Some(m_mle_proof), Some(r_mle_proof))
        } else {
            (None, None)
        };

        Proof {
            y,
            z1_sumcheck_proof,
            z3_sumcheck_proof,
            r_mle_proof,
            m_mle_proof,
        }
    }
}

fn compute_z1_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    m_mle: &MultilinearPolynomial<F>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    powers_of_alpha: &Vec<F>,
    vec_tau: &Vec<F>,
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
            ct.get_ring_element(1)
                .unwrap()
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
    vec_tau: &Vec<F>,
) -> Z3Mles<F> {
    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank

    let r_mle_base_prime_f_evals = vec_quotients
        .iter()
        .flat_map(|quotients| &quotients[1].coeffs[0..D]) // only consider the lower order coefficients
        .copied()
        .collect::<Vec<F::BasePrimeField>>();

    let r_mle_evals = r_mle_base_prime_f_evals
        .iter()
        .map(|e| F::from_base_prime_field(*e))
        .collect::<Vec<F>>();
    let r_mle_evals_len = r_mle_evals.len();

    let r_mle_over_base_prime_f = MultilinearPolynomial::new(r_mle_base_prime_f_evals, z3_num_vars);
    let r_mle = MultilinearPolynomial::new(r_mle_evals, z3_num_vars);

    let padded_alpha_mle_evals = (0..r_mle_evals_len / D)
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();
    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, z3_num_vars);

    // compute eq polynomial
    let m = (r_mle_evals_len / D).ilog2() as usize;
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << m);
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let padded_eq_tau_mle_evals = eq_tau_mle_evals
        .iter()
        .map(|e| vec![*e; D])
        .flatten()
        .collect::<Vec<F>>();

    let eq_tau_mle = MultilinearPolynomial::new(padded_eq_tau_mle_evals, z3_num_vars);

    Z3Mles {
        mles_over_f: vec![eq_tau_mle, r_mle, alpha_mle],
        r_mle_over_base_prime_f,
    }
}
