use std::marker::PhantomData;

use ark_ff::PrimeField;

use crate::{arith::linalg::Matrix, protocol::Proof, rlwe::RLWE};

pub struct Prover<const D: usize, F: PrimeField> {
    _pd: PhantomData<F>,
}

// TODO: support `x` as a vector of RLWE elements
impl<const D: usize, F: PrimeField> Prover<D, F> {
    pub fn prove(m: Matrix<F>, x: RLWE<D, F>) -> Proof<D, F> {
        let m_pr = m.process();

        // interpret each row as a ring elements
        let relts = m_pr.to_cyclotomic_relts::<D>();

        let y = relts
            .iter()
            .map(|pi| x.mul_constant(&pi))
            .collect::<Vec<RLWE<D, F>>>();

        Proof { y }
    }
}
