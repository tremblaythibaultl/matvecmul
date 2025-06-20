use std::marker::PhantomData;

use ark_ff::{CyclotomicMultSubgroup, PrimeField};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, ring::Ring},
    protocol::Proof,
    rlwe::RLWE,
};

pub struct Prover<const D: usize, F: PrimeField> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: PrimeField> Prover<D, F> {
    pub fn prove(m: &Matrix<F>, x: &Vec<RLWE<CyclotomicRing<D, F>>>) -> Proof<D, F> {
        // interpret each row as a ring elements
        let m_rq = m.lift_to_rq::<D>();

        println!("m_rq: {:#?}", m_rq);

        // mat vec mul
        let y = m_rq.mat_rlwe_vec_mul(&x);

        Proof { y }
    }
}
