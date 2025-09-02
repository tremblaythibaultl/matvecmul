use std::marker::PhantomData;

use ark_ff::{Field, PrimeField};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{Proof, sample_random_challenge, sumcheck::multilinear::MultilinearPolynomial},
    rlwe::{K, RLWE},
};

pub struct Prover<const D: usize, F: Field> {
    _pd: PhantomData<F>,
}

impl<const D: usize, F: Field> Prover<D, F> {
    pub fn prove(
        m: &Matrix<F::BasePrimeField>,
        x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    ) -> Proof<D, F::BasePrimeField> {
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

        let alpha = sample_random_challenge::<F>(true);
        let tau = sample_random_challenge::<F>(true);

        // println!("vec_remainders: {:#?}", vec_remainders);
        // println!("vec_quotients: {:#?}", vec_quotients);

        // (At least one of) the sumcheck polynomial(s) will be of degree 2 - so not multilinear.
        // I think it might be best to store several MLES and then perform the sumcheck for a product of MLES. Look at JolT implementation for more details.

        // this also clones the coefficients. should look into optimizing.
        let m_mle: MultilinearPolynomial<F::BasePrimeField> = m_rq.to_mle();

        // maybe there are better ways of computing this
        // or maybe we do not need to construct the polynomial and we can evaluate it on the fly
        let powers_of_alpha = (0..D).map(|i| alpha.pow([i as u64])).collect::<Vec<_>>();
        let alpha_mle = MultilinearPolynomial::new(
            powers_of_alpha.clone(),
            D.next_power_of_two().ilog2() as usize,
        );

        let x_alpha_mle_evals = x
            .iter()
            .map(|elem| {
                elem.get_ring_elements()
                    .into_iter()
                    .map(|poly| {
                        poly.coeffs
                            .iter()
                            .zip(powers_of_alpha.iter())
                            .map(|(c, a)| a.mul_by_base_prime_field(c))
                            .sum::<F>()
                    })
                    .collect::<Vec<F>>()
            })
            .flatten()
            .collect::<Vec<F>>();

        let x_alpha_mle: MultilinearPolynomial<F> = MultilinearPolynomial::new(
            x_alpha_mle_evals,
            (x.len() * (K + 1)).next_power_of_two().ilog2() as usize,
        );

        let r_alpha_mle_evals = vec_quotients
            .iter()
            .map(|elem| {
                elem.iter()
                    .map(|poly| {
                        poly.coeffs
                            .iter()
                            .zip(powers_of_alpha.iter())
                            .map(|(c, a)| a.mul_by_base_prime_field(c))
                            .sum::<F>()
                    })
                    .collect::<Vec<F>>()
            })
            .flatten()
            .collect::<Vec<F>>();

        let r_alpha_mle: MultilinearPolynomial<F> = MultilinearPolynomial::new(
            r_alpha_mle_evals,
            (vec_quotients.len() * vec_quotients[0].len())
                .next_power_of_two()
                .ilog2() as usize,
        );

        Proof {
            y,
            r: vec_quotients,
        }
    }
}
