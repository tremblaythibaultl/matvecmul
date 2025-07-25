use std::marker::PhantomData;

use ark_ff::PrimeField;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::Proof,
    rlwe::RLWE,
};

pub struct Prover<const D: usize, F: PrimeField> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: PrimeField> Prover<D, F> {
    pub fn prove(m: &Matrix<F>, x: &Vec<RLWE<CyclotomicRing<D, F>>>) -> Proof<D, F> {
        // interpret each row as a ring elements and rearrange columns
        let m_rq = m.lift_to_rq::<D>();

        // lift matrix entries to polynomial ring elements of max degree 2*D
        // this performs a clone of the coefficients
        // should look into optimizing
        let m_polyring = m_rq.lift_to_polynomial_ring();

        let x_polyring = x
            .iter()
            .map(|elem| elem.lift_to_polynomial_ring())
            .collect::<Vec<_>>();

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        //reduce each polynomial ring element to cyclotomic ring using long division
        let mut vec_quotients = Vec::<Vec<PolynomialRing<D, F>>>::with_capacity(y_polyring.len());
        let mut vec_remainders = Vec::<Vec<CyclotomicRing<D, F>>>::with_capacity(y_polyring.len());

        for elem in y_polyring {
            let (quotients, remainders) = elem.long_division_by_cyclotomic();
            vec_quotients.push(quotients);
            vec_remainders.push(remainders);
        }

        // mat vec mul -- this value should coincide with vec_remainders
        let y = m_rq.mat_rlwe_vec_mul(&x);

        println!("y: {:#?}", y);
        println!("vec_remainders: {:#?}", vec_remainders);

        Proof {
            y,
            r: vec_quotients,
        }
    }
}
