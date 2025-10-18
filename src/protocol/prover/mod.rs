use std::marker::PhantomData;

use ark_ff::{BigInteger, FftField, Field, PrimeField};
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

        // let poseidon_config = <F::BasePrimeField>::get_poseidon_config();
        // let mut transcript: PoseidonTranscript<F::BasePrimeField> =
        // PoseidonTranscript::new(poseidon_config);

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
    ) -> Proof<D, F> {
        let x_polyring = x
            .iter()
            .map(|elem| elem.lift_to_polynomial_ring())
            .collect::<Vec<_>>();

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        let (vec_quotients, _vec_remainders): (
            Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>,
            Vec<Vec<CyclotomicRing<D, F::BasePrimeField>>>,
        ) = y_polyring
            .par_iter()
            .map(|ct| ct.long_division_by_cyclotomic())
            .collect();

        // mat vec mul -- this value should coincide with vec_remainders
        let y = m_rq.mat_rlwe_vec_mul(&x);

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

        let y_bytes_to_absorb = y
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

        let bytes_to_absorb = &[x_bytes_to_absorb, y_bytes_to_absorb].concat();
        transcript.absorb_bytes(bytes_to_absorb);

        // TODO: absorb the commitment to r too

        let mut challenges = transcript.squeeze(m_rq.height().ilog2() as usize + 1);

        let alpha = challenges.pop().unwrap();
        let vec_tau = challenges;

        let powers_of_alpha: Vec<F> = (0..D)
            .scan(F::one(), |state, _| {
                let result = *state;
                *state *= alpha;
                Some(result)
            })
            .collect();

        let mut z1_mles =
            compute_z1_mles(&m_rq, &m_mle, x, z1_num_vars, &powers_of_alpha, &vec_tau);

        let z1_claim = sum_over_boolean_hypercube(&z1_mles);

        let (z1_sumcheck_proof, z1_challenges) = prove(z1_claim, &mut z1_mles, z1_num_vars);

        // We probably will be able to batch the two sumchecks (z_1 and z_3). Not clear how yet.

        // z_3
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let mut z3_mles = compute_z3_mles(&vec_quotients, &powers_of_alpha, z3_num_vars, &vec_tau);

        let z3_claim = sum_over_boolean_hypercube(&z3_mles.mles_over_f);

        let (z3_sumcheck_proof, z3_challenges) =
            prove(z3_claim, &mut z3_mles.mles_over_f, z3_num_vars);

        // Whir proofs for r_mle and m_mle.
        let mut rng = get_rng();
        let whir_r_mle = Whir::<F>::new(z3_mles.r_mle_over_base_prime_f.num_variables(), &mut rng);
        let r_mle_proof = whir_r_mle.prove(&z3_mles.r_mle_over_base_prime_f, &z3_challenges);

        let whir_m_mle = Whir::<F>::new(m_mle_over_base_f.num_variables(), &mut rng);
        let m_mle_proof = whir_m_mle.prove(&m_mle_over_base_f, &z1_challenges);

        println!("m_mle_proof_size: {}", m_mle_proof.proof.len());

        // sanity check with z_2

        // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank
        // let y_alpha_mle_evals = y
        //     .iter()
        //     .map(|ct| {
        //         ct.get_ring_elements()[1]
        //             .coeffs
        //             .iter()
        //             .zip(powers_of_alpha.iter())
        //             .map(|(c, a)| a.mul_by_base_prime_field(c))
        //             .sum::<F>()
        //     })
        //     .collect::<Vec<F>>();

        // let z2_num_vars = m_rq.height().ilog2() as usize;
        // let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << z2_num_vars);
        // build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

        // let _z_2 = sum_over_boolean_hypercube(&[
        //     MultilinearPolynomial::new(y_alpha_mle_evals, z2_num_vars),
        //     MultilinearPolynomial::new(eq_tau_mle_evals, z2_num_vars),
        // ]);

        Proof {
            y,
            r: vec_quotients,
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
