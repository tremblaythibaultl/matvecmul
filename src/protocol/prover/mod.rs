use std::marker::PhantomData;

use ark_ff::{FftField, Field};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

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
            // consider body only (compressed ciphertext)
            let el = ct.get_ring_element(1).unwrap();
            el.serialize(&mut x_bytes_to_absorb).unwrap();
        }

        transcript.absorb_bytes_par(&x_bytes_to_absorb);

        let mut beta = transcript.squeeze(1);
        beta.push(beta[0] * beta[0]);

        let x_polyring = x
            .iter()
            .map(|elem| elem.lift_to_polynomial_ring())
            .collect::<Vec<_>>();

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        let (vec_quotients, _vec_remainders): (
            Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>,
            Vec<Vec<CyclotomicRing<D, F::BasePrimeField>>>,
        ) = y_polyring
            .iter()
            .map(|ct| ct.long_division_by_cyclotomic())
            .collect();

        let mut rng = get_rng();
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let (r0_mle_base_prime_f_evals, r1_mle_base_prime_f_evals) =
            compute_r_bpf_mle_evals::<D, F>(&vec_quotients);

        let r0_mle_over_base_prime_f =
            MultilinearPolynomial::new(r0_mle_base_prime_f_evals.clone(), z3_num_vars);
        let whir_r0_mle = Whir::<F>::new(r0_mle_over_base_prime_f.num_variables(), &mut rng);

        let r1_mle_over_base_prime_f =
            MultilinearPolynomial::new(r1_mle_base_prime_f_evals.clone(), z3_num_vars);
        let whir_r1_mle = Whir::<F>::new(r1_mle_over_base_prime_f.num_variables(), &mut rng);

        let mut r0_commitment = vec![];
        let whir_r0_mle_commitment_and_prover_state = if include_pcs {
            let (commitment, prover_state) = whir_r0_mle.commit(&r0_mle_over_base_prime_f);
            r0_commitment.extend_from_slice(prover_state.narg_string());
            Some((commitment, prover_state))
        } else {
            None
        };

        let mut r1_commitment = vec![];
        let whir_r1_mle_commitment_and_prover_state = if include_pcs {
            let (commitment, prover_state) = whir_r1_mle.commit(&r1_mle_over_base_prime_f);
            r1_commitment.extend_from_slice(prover_state.narg_string());
            Some((commitment, prover_state))
        } else {
            None
        };

        let whir_m_mle = Whir::<F>::new(m_mle_over_base_f.num_variables(), &mut rng);

        let mut m_commitment = vec![];
        let whir_m_mle_commitment_and_prover_state = if include_pcs {
            let (commitment, prover_state) = whir_m_mle.commit(&m_mle_over_base_f);
            m_commitment.extend_from_slice(prover_state.narg_string());
            Some((commitment, prover_state))
        } else {
            None
        };

        // mat vec mul -- this value should coincide with vec_remainders
        // maybe it would be better to compute it from `y_polyring`
        let y = m_rq.mat_rlwe_vec_mul(&x);

        // assumes RLWE rank is 1
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

        let mut y_bytes_to_absorb = vec![];
        for y_batched_i in y_batched.iter() {
            y_batched_i
                .serialize_uncompressed(&mut y_bytes_to_absorb)
                .unwrap();
        }

        transcript.absorb_bytes_par(&y_bytes_to_absorb);
        transcript.absorb_bytes(&r0_commitment);
        transcript.absorb_bytes(&r1_commitment);

        let mut challenges = transcript.squeeze(m_rq.height().ilog2() as usize + 1);

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

        let mut z1_mles = compute_z1_mles(
            &m_rq,
            &m_mle,
            x,
            z1_num_vars,
            &beta,
            &powers_of_alpha,
            &vec_tau,
        );

        let z1_claim = sum_over_boolean_hypercube(&z1_mles);

        let (z1_sumcheck_proof, z1_challenges) =
            prove(z1_claim, &mut z1_mles, z1_num_vars, transcript);

        // z_3
        let mut z3_mles = compute_z3_mles::<D, F>(
            &powers_of_alpha,
            &r0_mle_base_prime_f_evals,
            &r1_mle_base_prime_f_evals,
            &beta,
            z3_num_vars,
            &vec_tau,
        );

        let z3_claim = sum_over_boolean_hypercube(&z3_mles);

        let (z3_sumcheck_proof, z3_challenges) =
            prove(z3_claim, &mut z3_mles, z3_num_vars, transcript);

        // Whir proofs for r_mle and m_mle.
        let (m_mle_proof, r0_mle_proof, r1_mle_proof) = if include_pcs {
            let (r0_mle_commitment, r0_mle_prover_state) =
                whir_r0_mle_commitment_and_prover_state.unwrap();

            let r0_mle_proof =
                whir_r0_mle.prove(r0_mle_commitment, r0_mle_prover_state, &z3_challenges);

            let (r1_mle_commitment, r1_mle_prover_state) =
                whir_r1_mle_commitment_and_prover_state.unwrap();

            let r1_mle_proof =
                whir_r1_mle.prove(r1_mle_commitment, r1_mle_prover_state, &z3_challenges);

            let (m_mle_commitment, m_mle_prover_state) =
                whir_m_mle_commitment_and_prover_state.unwrap();
            let m_mle_proof =
                whir_m_mle.prove(m_mle_commitment, m_mle_prover_state, &z1_challenges);

            (Some(m_mle_proof), Some(r0_mle_proof), Some(r1_mle_proof))
        } else {
            (None, None, None)
        };

        Proof {
            y,
            z1_sumcheck_proof,
            z3_sumcheck_proof,
            r0_commitment,
            r1_commitment,
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
