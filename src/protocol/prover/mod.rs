use std::marker::PhantomData;

use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        Proof, compute_mles, preprocess, sample_random_challenge,
        sumcheck::{multilinear::MultilinearPolynomial, prove},
        utils::sum_over_boolean_hypercube,
    },
    rlwe::RLWE,
};

pub struct Prover<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: Field> Prover<D, F> {
    pub fn prove(
        m: &Matrix<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    ) -> Proof<D, F> {
        let (m_rq, m_polyring, x_polyring, m_mle, num_vars) = preprocess::<D, F>(m, x);

        let y_polyring = m_polyring.mat_rlwe_vec_mul(&x_polyring);

        //reduce each polynomial ring element to cyclotomic ring using long division
        let mut vec_quotients =
            Vec::<Vec<PolynomialRing<D, F::BasePrimeField>>>::with_capacity(y_polyring.len());
        let mut vec_remainders =
            Vec::<Vec<CyclotomicRing<D, F::BasePrimeField>>>::with_capacity(y_polyring.len());

        for elem in y_polyring {
            let (quotients, remainders) = elem.long_division_by_cyclotomic();
            vec_quotients.push(quotients);
            vec_remainders.push(remainders);
        }

        // mat vec mul -- this value should coincide with vec_remainders
        let y = m_rq.mat_rlwe_vec_mul(&x);

        println!("y: {:#?}", y);
        println!("vec_remainders: {:#?}", vec_remainders);
        println!("vec_quotients: {:#?}", vec_quotients);

        // TODO: implement fiat shamir
        let alpha = sample_random_challenge::<F>(true);
        let tau = sample_random_challenge::<F>(true);

        let mut mles = compute_mles(&m_rq, &m_mle, x, num_vars, &alpha, &tau);

        let claim = sum_over_boolean_hypercube(&mles);

        let (z1_sumcheck_proof, challenges) = prove(claim, &mut mles, num_vars);

        // We probably will be able to batch the two sumchecks (z_1 and z_3). Not clear how yet.

        // z_3
        // let r_alpha_mle_evals = vec_quotients
        //     .iter()
        //     .map(|elem| {
        //         elem.iter()
        //             .map(|poly| {
        //                 poly.coeffs
        //                     .iter()
        //                     .zip(powers_of_alpha.iter())
        //                     .map(|(c, a)| a.mul_by_base_prime_field(c))
        //                     .sum::<F>()
        //             })
        //             .collect::<Vec<F>>()
        //     })
        //     .flatten()
        //     .collect::<Vec<F>>();

        // let r_alpha_mle: MultilinearPolynomial<F> = MultilinearPolynomial::new(
        //     r_alpha_mle_evals,
        //     (vec_quotients.len() * vec_quotients[0].len())
        //         .next_power_of_two()
        //         .ilog2() as usize,
        // );

        Proof {
            y,
            r: vec_quotients,
            z1_sumcheck_proof,
        }
    }
}
