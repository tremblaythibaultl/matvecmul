use ark_ff::PrimeField;

use crate::{arith::cyclotomic_ring::CyclotomicRing, rlwe::RLWE};

pub mod prover;
mod verifier;

pub struct Proof<const D: usize, F: PrimeField> {
    pub y: Vec<RLWE<CyclotomicRing<D, F>>>,
}
