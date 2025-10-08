use std::marker::PhantomData;

use ark_crypto_primitives::{
    crh::{CRHScheme, TwoToOneCRHScheme},
    merkle_tree::Config,
};
use ark_ff::{FftField, Field};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::RngCore;
use spongefish::{DomainSeparator, ProverState, VerifierState};
use spongefish_pow::blake3::Blake3PoW;
use whir::{
    crypto::merkle_tree::{
        blake3::{Blake3Compress, Blake3LeafHash, Blake3MerkleTreeParams},
        parameters::default_config,
    },
    parameters::{
        DeduplicationStrategy, FoldingFactor, MerkleProofStrategy, MultivariateParameters,
        ProtocolParameters, SoundnessType, default_max_pow,
    },
    poly_utils::{evals::EvaluationsList, multilinear::MultilinearPoint},
    whir::{
        committer::{CommitmentReader, CommitmentWriter},
        domainsep::{DigestDomainSeparator, WhirDomainSeparator},
        parameters::WhirConfig,
        prover::Prover,
        statement::{Statement, Weights},
        utils::{DigestToUnitDeserialize, DigestToUnitSerialize},
        verifier::Verifier,
    },
};

use crate::protocol::sumcheck::multilinear::MultilinearPolynomial;

pub type PowStrategy = Blake3PoW;

#[derive(Clone)]
pub struct WhirProof<F> {
    pub claim: F,
    pub proof: Vec<u8>,
}

pub struct Whir<
    F: FftField + CanonicalSerialize,
    MerkleConfig: Config<Leaf = [F]> + Clone = Blake3MerkleTreeParams<F>,
> {
    num_variables: usize,
    domainsep: DomainSeparator,
    params: WhirConfig<F, MerkleConfig, PowStrategy>,
    _phantom: PhantomData<MerkleConfig>,
}

impl<F, MerkleConfig> Whir<F, MerkleConfig>
where
    F: Field + FftField + CanonicalSerialize,
    MerkleConfig: Config<Leaf = [F]> + Clone,
    MerkleConfig::InnerDigest: AsRef<[u8]> + From<[u8; 32]>,
    DomainSeparator: DigestDomainSeparator<MerkleConfig>,
    ProverState: DigestToUnitSerialize<MerkleConfig>,
    for<'a> VerifierState<'a>: DigestToUnitDeserialize<MerkleConfig>,
{
    // TODO: Parameters need to be revised for performance/security.
    pub const SECURITY_LEVEL: usize = 80;
    pub const RATE: usize = 1;
    pub const FIRST_ROUND_FOLDING_FACTOR: usize = 4;
    pub const FOLDING_FACTOR: usize = 4;
    pub const SOUNDNESS_TYPE: SoundnessType = SoundnessType::UniqueDecoding;

    pub fn new<R: RngCore>(num_variables: usize, rng: &mut R) -> Self
    where
        <<MerkleConfig as Config>::LeafHash as CRHScheme>::Parameters: From<()>,
        <<MerkleConfig as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters: From<()>,
    {
        let (leaf_hash_params, two_to_one_params) =
            default_config::<F, Blake3LeafHash<F>, Blake3Compress>(rng);

        let multivariate_params = MultivariateParameters::<F>::new(num_variables);
        let whir_params = ProtocolParameters::<MerkleConfig, PowStrategy> {
            initial_statement: true,
            security_level: Self::SECURITY_LEVEL,
            pow_bits: default_max_pow(num_variables, Self::RATE),
            folding_factor: FoldingFactor::ConstantFromSecondRound(
                Self::FIRST_ROUND_FOLDING_FACTOR,
                Self::FOLDING_FACTOR,
            ),
            leaf_hash_params: leaf_hash_params.into(),
            two_to_one_params: two_to_one_params.into(),
            soundness_type: Self::SOUNDNESS_TYPE,
            _pow_parameters: Default::default(),
            starting_log_inv_rate: Self::RATE,
            batch_size: 1,
            deduplication_strategy: DeduplicationStrategy::Enabled,
            merkle_proof_strategy: MerkleProofStrategy::Compressed,
        };

        let params =
            WhirConfig::<F, MerkleConfig, PowStrategy>::new(multivariate_params, whir_params);

        Self {
            num_variables,
            domainsep: DomainSeparator::new("matvecmul")
                .commit_statement(&params)
                .add_whir_proof(&params),
            params,
            _phantom: PhantomData,
        }
    }

    // Proves that `poly` evaluates to a particular value at `point`.
    pub fn prove(&self, poly: &MultilinearPolynomial<F>, point: &[F]) -> WhirProof<F> {
        let mut prover_state = self.domainsep.to_prover_state();
        // We assume that the polynomial is defined over the base prime field (e.g. for r_mle in compute_z3_mles()).
        // Otherwise, we could define a MultilinearPolynomial<F::BasePrimeField> type for this purpose and refactor code a bit.
        let base_prime_evals: Vec<F::BasePrimeField> = poly
            .evals()
            .iter()
            .map(|e| {
                let mut components = e.to_base_prime_field_elements();
                let first = components.next().unwrap();
                assert_eq!(components.next(), None);
                first
            })
            .collect();
        let evals = EvaluationsList::new(base_prime_evals);
        let coeffs = evals.to_coeffs();
        let committer = CommitmentWriter::new(self.params.clone());
        let witness = committer.commit(&mut prover_state, &coeffs).unwrap();
        let mut statement: Statement<F> = Statement::<F>::new(self.num_variables);
        let point = MultilinearPoint(point.to_vec());
        let claim = coeffs.evaluate_at_extension(&point);
        let weights = Weights::evaluation(point.clone());
        statement.add_constraint(weights, claim);
        let prover = Prover::new(self.params.clone());

        prover
            .prove(&mut prover_state, statement.clone(), witness)
            .unwrap();

        WhirProof {
            proof: prover_state.narg_string().to_vec(),
            claim,
        }
    }

    // Returns Ok if proof verification succeeded and Err otherwise.
    pub fn verify(&self, proof: &WhirProof<F>, point: &[F]) -> anyhow::Result<()> {
        let commitment_reader = CommitmentReader::new(&self.params);
        let verifier = Verifier::new(&self.params);
        let mut verifier_state = self.domainsep.to_verifier_state(&proof.proof);
        let parsed_commitment = commitment_reader.parse_commitment(&mut verifier_state)?;
        let mut statement: Statement<F> = Statement::<F>::new(self.num_variables);
        let point = MultilinearPoint(point.to_vec());
        let weights = Weights::evaluation(point.clone());
        statement.add_constraint(weights, proof.claim);
        verifier.verify(&mut verifier_state, &parsed_commitment, &statement)?;
        Ok(())
    }
}
