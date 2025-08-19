use std::marker::PhantomData;

use ark_ff::Field;

use crate::protocol::Proof;

pub struct Verifier<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: Field> Verifier<D, F> {
    fn verify(pr: Proof<D, F>) -> bool {
        todo!()
    }
}
