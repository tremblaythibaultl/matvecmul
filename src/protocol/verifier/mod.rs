use std::marker::PhantomData;

use ark_ff::{FftField, Field};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, field::GetPoseidonConfig, linalg::Matrix},
    protocol::{
        Proof, pcs::whir::Whir, sumcheck::multilinear::MultilinearPolynomial,
        transcript::PoseidonTranscript, utils::build_eq_poly,
    },
    rand::get_rng,
    rlwe::RLWE,
};

pub struct Verifier<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F> Verifier<D, F>
where
    F: FftField,
    F::BasePrimeField: GetPoseidonConfig<F::BasePrimeField> + ark_crypto_primitives::sponge::Absorb,
{
    pub fn preprocess(
        m: &Matrix<F::BasePrimeField>,
    ) -> (
        Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        MultilinearPolynomial<F>,
        usize,
        PoseidonTranscript<F::BasePrimeField>,
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

        let poseidon_config = <F::BasePrimeField>::get_poseidon_config();
        let mut transcript: PoseidonTranscript<F::BasePrimeField> =
            PoseidonTranscript::new(poseidon_config);

        for elem in m.data.iter() {
            transcript.absorb(elem);
        }

        (m_rq, m_mle, num_vars, transcript)
    }

    pub fn verify(
        m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
        m_mle: &MultilinearPolynomial<F>,
        z1_num_vars: usize,
        transcript: &mut PoseidonTranscript<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
        proof: Proof<D, F>,
    ) -> Result<Vec<F>, ()> {
        // only consider mask for now
        for ct in x.iter() {
            for coeff in ct.get_ring_elements()[0].coeffs.iter() {
                transcript.absorb(coeff);
            }
        }

        // only consider mask for now
        for ct in proof.y.iter() {
            for coeff in ct.get_ring_elements()[0].coeffs.iter() {
                transcript.absorb(coeff);
            }
        }

        let alpha =
            F::from_base_prime_field_elems([transcript.squeeze(), transcript.squeeze()]).unwrap();

        let m = m_rq.height().ilog2() as usize;
        let mut vec_tau = Vec::<F>::with_capacity(m);
        for _ in 0..m {
            vec_tau.push(
                F::from_base_prime_field_elems([transcript.squeeze(), transcript.squeeze()])
                    .unwrap(),
            );
        }

        // compute unpadded MLEs
        let unpadded_z1_mles = compute_z1_mles(&m_rq, x, z1_num_vars, &alpha, &vec_tau);

        let powers_of_alpha = unpadded_z1_mles[2].evals();

        let (z1_original_claim, z1_final_claim, z1_challenges) =
            proof.z1_sumcheck_proof.verify(z1_num_vars, 4).unwrap();

        let eq_eval = evaluate_eq_tau_at_challenges(&unpadded_z1_mles[0], &z1_challenges[0..m]);

        let x_alpha_eval = evaluate_x_alpha_at_challenges(
            &unpadded_z1_mles[1],
            &z1_challenges[m_rq.height().ilog2() as usize
                ..m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize]
                .to_vec(),
        );

        let z1_ell_eval = evaluate_ell_at_challenges(
            &unpadded_z1_mles[2],
            &z1_challenges[m_rq.height().ilog2() as usize + m_rq.width().ilog2() as usize..],
        );

        // This does NOT need to be computed by the verifier
        // The value should come from the PCS.
        // Included for test purposes
        let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << z1_num_vars);
        build_eq_poly(&z1_challenges, &mut eq_x_challenges_mle_evals);
        let m_eval = m_mle
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        let z1_eval_at_random_point = eq_eval * m_eval * x_alpha_eval * z1_ell_eval;

        if z1_eval_at_random_point != z1_final_claim {
            return Err(());
        } else {
            println!("z1 Verification successful");
        }

        // z_3 verification
        let z3_num_vars = z1_num_vars - m_rq.width().ilog2() as usize;

        let (z3_original_claim, z3_final_claim, z3_challenges) =
            proof.z3_sumcheck_proof.verify(z3_num_vars, 3).unwrap();

        let z3_ell_eval = evaluate_ell_at_challenges(
            &unpadded_z1_mles[2],
            &z3_challenges[m_rq.height().ilog2() as usize..].to_vec(),
        );

        // This does NOT need to be computed by the verifier
        // The value should come from the PCS.
        // Included for test purposes
        let mut eq_x_challenges_mle_evals = Vec::<F>::with_capacity(1 << z3_num_vars);
        build_eq_poly(&z3_challenges, &mut eq_x_challenges_mle_evals);
        let r_eval = proof
            .r_mle
            .evals()
            .iter()
            .zip(eq_x_challenges_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        let mut rng = get_rng();
        let whir = Whir::<F>::new(z3_num_vars, &mut rng);
        whir.verify(&proof.r_mle_proof, &z3_challenges)
            .expect("r_mle proof does not verify");
        assert_eq!(r_eval, proof.r_mle_proof.claim);

        // eq_tau eval
        let eq_eval = evaluate_eq_tau_at_challenges(&unpadded_z1_mles[0], &z1_challenges[0..m]);

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

        let mut eq_tau_mle_evals = Vec::<F>::with_capacity(1 << z2_num_vars);
        build_eq_poly(&vec_tau[0..z2_num_vars], &mut eq_tau_mle_evals);

        // focus on the first component of the ciphertext only for now.
        let y_alpha_mle_evals = proof
            .y
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

        assert_eq!(y_alpha_mle_evals.len(), eq_tau_mle_evals.len());

        let z_2 = y_alpha_mle_evals
            .iter()
            .zip(eq_tau_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

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
