use ark_ff::PrimeField;

use crate::rlwe::RLWE;

pub mod prover;
mod verifier;

pub struct Proof<const D: usize, F: PrimeField> {
    pub y: Vec<RLWE<D, F>>,
}
