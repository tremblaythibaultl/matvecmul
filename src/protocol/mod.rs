use ark_ff::Field;

use crate::{
    arith::cyclotomic_ring::CyclotomicRing,
    protocol::{pcs::whir::WhirProof, sumcheck::SumCheckProof},
    rlwe::RLWE,
};

pub mod pcs;
pub mod prover;
pub mod sumcheck;
pub mod transcript;
mod utils;
pub mod verifier;

#[derive(Clone)]
pub struct Proof<const D: usize, F: Field> {
    pub y: Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    pub z1_sumcheck_proof: SumCheckProof<F>,
    pub z3_sumcheck_proof: SumCheckProof<F>,
    pub r0_commitment: Vec<u8>,
    pub r1_commitment: Vec<u8>,
    pub r0_mle_proof: Option<WhirProof<F>>,
    pub r1_mle_proof: Option<WhirProof<F>>,
    pub m_mle_proof: Option<WhirProof<F>>,
}
