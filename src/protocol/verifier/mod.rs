use std::marker::PhantomData;

use ark_ff::PrimeField;

use crate::protocol::Proof;

pub struct Verifier<const D: usize, F: PrimeField> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: PrimeField> Verifier<D, F> {
    fn verify(pr: Proof<D, F>) -> bool {
        todo!()
    }
}
