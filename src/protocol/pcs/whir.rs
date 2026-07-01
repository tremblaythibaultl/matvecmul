use ark_ff::Field;
use ark_std::rand::{
    RngCore,
    distributions::{Distribution, Standard},
};
use std::borrow::Cow;
use whir::{
    algebra::{
        embedding::Basefield,
        linear_form::{Evaluate, LinearForm, MultilinearExtension},
    },
    hash,
    parameters::ProtocolParameters,
    protocols::whir::Config,
    transcript::{Codec, DomainSeparator, Proof, ProverState, VerifierState, codecs::Empty},
};

use crate::protocol::sumcheck::multilinear::MultilinearPolynomial;

#[derive(Clone)]
pub struct WhirProof<F: Field> {
    pub claim: F,
    pub proof: Proof,
}

impl<F: Field> WhirProof<F> {
    pub fn size_in_bytes(&self) -> usize {
        self.proof.narg_string.len() + self.proof.hints.len()
    }
}

pub struct WhirCommitment<F: Field> {
    evals: Vec<F::BasePrimeField>,
    commitment_bytes: Vec<u8>,
}

impl<F: Field> WhirCommitment<F> {
    pub fn narg_string(&self) -> &[u8] {
        &self.commitment_bytes
    }
}

pub struct Whir<F: Field> {
    num_variables: usize,
    config: Config<Basefield<F>>,
}

impl<F> Whir<F>
where
    F: Field,
    Standard: Distribution<F> + Distribution<F::BasePrimeField>,
    F: Codec<[u8]>,
{
    pub const SECURITY_LEVEL: usize = 100;
    pub const RATE: usize = 1;
    pub const FIRST_ROUND_FOLDING_FACTOR: usize = 4;
    pub const FOLDING_FACTOR: usize = 4;
    pub const BATCH_SIZE: usize = 1;

    pub fn new<R: RngCore>(num_variables: usize, _rng: &mut R) -> Self {
        let whir_params = ProtocolParameters {
            security_level: Self::SECURITY_LEVEL,
            // Mirrors the previous `default_max_pow(num_variables, rate)`,
            // which was `num_variables + log_inv_rate - 3`.
            pow_bits: (num_variables + Self::RATE).saturating_sub(3),
            initial_folding_factor: Self::FIRST_ROUND_FOLDING_FACTOR,
            folding_factor: Self::FOLDING_FACTOR,
            unique_decoding: false,
            starting_log_inv_rate: Self::RATE,
            batch_size: Self::BATCH_SIZE,
            hash_id: hash::BLAKE3,
        };

        let config = Config::<Basefield<F>>::new(1 << num_variables, &whir_params);

        Self {
            num_variables,
            config,
        }
    }

    fn whir_point(point: &[F]) -> Vec<F> {
        point.to_vec()
    }

    fn domain_separator(&self) -> DomainSeparator<'static, Empty> {
        DomainSeparator::protocol(&self.config)
            .session(&"matvecmul")
            .instance(&Empty)
    }

    pub fn commit(&self, poly: &MultilinearPolynomial<F::BasePrimeField>) -> WhirCommitment<F> {
        let evals = poly.evals().to_vec();

        let ds = self.domain_separator();
        let mut prover_state = ProverState::new_std(&ds);
        let _witness = self.config.commit(&mut prover_state, &[&evals]);
        let commitment_bytes = prover_state.proof().narg_string;

        WhirCommitment {
            evals,
            commitment_bytes,
        }
    }

    pub fn prove(&self, commitment: WhirCommitment<F>, point: &[F]) -> WhirProof<F> {
        let WhirCommitment { evals, .. } = commitment;

        let ds = self.domain_separator();
        let mut prover_state = ProverState::new_std(&ds);
        let witness = self.config.commit(&mut prover_state, &[&evals]);

        let form = MultilinearExtension::new(Self::whir_point(point));
        let claim = form.evaluate(self.config.embedding(), &evals);

        let _ = self.config.prove(
            &mut prover_state,
            vec![Cow::Borrowed(evals.as_slice())],
            vec![Cow::Owned(witness)],
            vec![Box::new(form) as Box<dyn LinearForm<F>>],
            Cow::Owned(vec![claim]),
        );

        WhirProof {
            claim,
            proof: prover_state.proof(),
        }
    }

    // Proves that `poly` evaluates to a particular value at `point`.
    pub fn commit_and_prove(
        &self,
        poly: &MultilinearPolynomial<F::BasePrimeField>,
        point: &[F],
    ) -> WhirProof<F> {
        let evals = poly.evals().to_vec();

        let ds = self.domain_separator();
        let mut prover_state = ProverState::new_std(&ds);
        let witness = self.config.commit(&mut prover_state, &[&evals]);

        let form = MultilinearExtension::new(Self::whir_point(point));
        let claim = form.evaluate(self.config.embedding(), &evals);

        let _ = self.config.prove(
            &mut prover_state,
            vec![Cow::Borrowed(evals.as_slice())],
            vec![Cow::Owned(witness)],
            vec![Box::new(form) as Box<dyn LinearForm<F>>],
            Cow::Owned(vec![claim]),
        );

        WhirProof {
            claim,
            proof: prover_state.proof(),
        }
    }

    // Returns Ok if proof verification succeeded and Err otherwise.
    pub fn verify(&self, proof: &WhirProof<F>, point: &[F]) -> anyhow::Result<()> {
        let _ = self.num_variables;
        let ds = self.domain_separator();
        let mut verifier_state = VerifierState::new_std(&ds, &proof.proof);

        let commitment = self
            .config
            .receive_commitment(&mut verifier_state)
            .map_err(|_| anyhow::anyhow!("failed to parse WHIR commitment"))?;

        let evaluations = [proof.claim];
        let final_claim = self
            .config
            .verify(&mut verifier_state, &[&commitment], &evaluations)
            .map_err(|_| anyhow::anyhow!("WHIR proof verification failed"))?;

        // Tie the deferred multilinear-extension constraint back to the claim.
        let form = MultilinearExtension::new(Self::whir_point(point));
        final_claim
            .verify([&form as &dyn LinearForm<F>])
            .map_err(|_| {
                anyhow::anyhow!("WHIR final claim does not match the committed polynomial")
            })?;

        Ok(())
    }

    pub fn num_variables(&self) -> usize {
        self.num_variables
    }
}
