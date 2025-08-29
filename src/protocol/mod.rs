use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, polynomial_ring::PolynomialRing},
    rlwe::RLWE,
};

pub mod prover;
pub mod sumcheck;
mod verifier;

pub struct Proof<const D: usize, F: Field> {
    pub y: Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    pub r: Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>, // should be a commitment to `r`
}

pub fn sample_random_challenge<F: Field>(test_mode: bool) -> F {
    if test_mode {
        // In test mode, return a fixed challenge

        // let coeffs = (1..=F::extension_degree())
        //     .map(|i| F::BasePrimeField::from(i as u64))
        //     .collect::<Vec<_>>();

        // return F::from_base_prime_field_elems(coeffs).unwrap(); // Should be safe to unwrap as the slice length matches the extension degree
        return F::from(2);
    } else {
        // TODO: implement a proper random challenge sampling
        // Should basically hash everything that is known to the verifier, i.e. `M`, `x`, `y`, and `r`.
        // The random challenge will probably be an element of `Gal(p, d)` for some `d` (probably 2?) so we should choose an appropriate hash function (Poseidon?)
        todo!()
    }
}
