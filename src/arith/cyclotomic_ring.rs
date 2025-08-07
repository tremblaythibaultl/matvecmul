use std::{
    iter::Sum,
    sync::{Once, OnceLock},
};

use ark_ff::{Field, PrimeField};
use ark_std::rand::Rng;

use tfhe_ntt::prime64::Plan;

use super::ring::Ring;

// At what point we start using NTT for multiplication.
const NTT_CUTOFF: usize = 32;

/// Represents an element of F\[X\]/(X^N + 1).
#[derive(Clone, Debug)]
pub struct CyclotomicRing<const D: usize, F: PrimeField> {
    pub coeffs: Vec<F>,
}

impl<const D: usize, F: PrimeField> CyclotomicRing<D, F> {
    pub fn from_coeffs(coeffs: &[F]) -> Self {
        assert!(
            coeffs.len() == D,
            "Coefficient vector must have length equal to the degree of the cyclotomic ring"
        );

        CyclotomicRing {
            coeffs: coeffs.to_vec(),
        }
    }

    pub fn add_constant(&self, constant: F) -> Self {
        let mut res: CyclotomicRing<D, F> = self.clone();
        res.coeffs[0] = res.coeffs[0] + constant;
        res
    }

    pub fn add_constant_assign(&mut self, constant: F) {
        self.coeffs[0] = self.coeffs[0] + constant;
    }

    /// Multiplies the residue polynomial by X^{exponent} = X^{2N + exponent}.
    /// `exponent` is assumed to be reduced modulo 2N.
    pub fn multiply_by_monomial(&self, exponent: usize) -> Self {
        let mut rotated_coeffs = Vec::<F>::with_capacity(D);

        let reverse = exponent >= D;
        let exponent = exponent % D;

        for i in 0..D {
            rotated_coeffs.push({
                if i < exponent {
                    if reverse {
                        self.coeffs[i + D - exponent]
                    } else {
                        -self.coeffs[i + D - exponent]
                    }
                } else if reverse {
                    -self.coeffs[i - exponent]
                } else {
                    self.coeffs[i - exponent]
                }
            })
        }

        CyclotomicRing {
            coeffs: rotated_coeffs,
        }
    }

    pub fn get_random() -> Self {
        let mut rng = ark_std::test_rng();

        let coeffs = (0..D).map(|_| F::rand(&mut rng)).collect();

        Self { coeffs }
    }

    pub fn get_random_bin() -> Self {
        let mut rng = ark_std::test_rng();

        let coeffs = (0..D).map(|_| F::from(rng.gen_range(0..=1))).collect();

        Self { coeffs }
    }

    pub fn basic_mul(&self, rhs: &Self) -> Self {
        let mut coeffs = Vec::<F>::with_capacity(D);
        for i in 0..D {
            let mut coeff = F::ZERO;
            for j in 0..i + 1 {
                coeff = coeff + (self.coeffs[j] * rhs.coeffs[i - j]);
            }
            for j in i + 1..D {
                coeff = coeff - (self.coeffs[j] * rhs.coeffs[D - j + i]);
            }
            coeffs.push(coeff);
        }
        CyclotomicRing { coeffs }
    }

    // If the plan is valid for 64 bit NTT based on D and F, we use NTT. Else we use basic multiplication.
    fn try_ntt_mul(&self, rhs: &Self) -> Self {
        if let Some(plan) = Self::ntt_plan_64() {
            self.ntt_mul(rhs, plan)
        } else {
            self.basic_mul(rhs)
        }
    }

    fn ntt_mul(&self, rhs: &Self, plan: &Plan) -> Self {
        // From https://github.com/zama-ai/tfhe-rs/blob/main/tfhe-ntt/examples/mul_poly_prime.rs
        // Converting to u64 should be safe as we assume that the modulus fits in 64 bits and that's checked at NTT plan creation.
        let mut lhs_u64 = self
            .coeffs
            .iter()
            .map(|c| c.into_bigint().as_ref()[0])
            .collect::<Vec<_>>();
        let mut rhs_u64 = rhs
            .coeffs
            .iter()
            .map(|c| c.into_bigint().as_ref()[0])
            .collect::<Vec<_>>();
        plan.fwd(&mut lhs_u64);
        plan.fwd(&mut rhs_u64);
        plan.mul_assign_normalize(&mut lhs_u64, &rhs_u64);
        plan.inv(&mut lhs_u64);
        let coeffs = lhs_u64.iter().map(|&x| F::from(x)).collect::<Vec<_>>();
        CyclotomicRing { coeffs }
    }

    // Gets a static plan for the 64 bit NTT. If the polynomial size or the prime modulus are not suitable,
    // returns None.
    fn ntt_plan_64() -> Option<&'static Plan> {
        static PLAN: OnceLock<Option<Plan>> = OnceLock::new();
        PLAN.get_or_init(|| {
            if F::MODULUS.as_ref().len() != 1 {
                // Modulus must fit in 64 bits.
                return None;
            }
            let modulus = F::MODULUS.as_ref()[0];
            Plan::try_new(D, modulus)
        })
        .as_ref()
    }
}

impl<const D: usize, F: PrimeField> Default for CyclotomicRing<D, F> {
    fn default() -> Self {
        CyclotomicRing {
            coeffs: vec![F::ZERO; D],
        }
    }
}

impl<const D: usize, F: PrimeField> Sum for CyclotomicRing<D, F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(CyclotomicRing::default(), |sum, summand| sum.add(&summand))
    }
}

impl<const D: usize, F: PrimeField> Ring for CyclotomicRing<D, F> {
    const DEGREE: usize = D;
    type BaseField = F;

    fn add(&self, rhs: &CyclotomicRing<D, F>) -> Self {
        let mut res = Self::default();
        for i in 0..D {
            res.coeffs[i] = self.coeffs[i] + rhs.coeffs[i];
        }
        res
    }

    fn add_assign(&mut self, rhs: &CyclotomicRing<D, F>) {
        for i in 0..D {
            self.coeffs[i] = self.coeffs[i] + rhs.coeffs[i];
        }
    }

    fn sub(&self, rhs: &CyclotomicRing<D, F>) -> Self {
        let mut res = Self::default();
        for i in 0..D {
            res.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
        res
    }

    fn sub_assign(&mut self, rhs: &CyclotomicRing<D, F>) {
        for i in 0..D {
            self.coeffs[i] = self.coeffs[i] - rhs.coeffs[i];
        }
    }

    fn mul(&self, rhs: &Self) -> Self {
        if D <= NTT_CUTOFF {
            return self.basic_mul(rhs);
        }
        self.try_ntt_mul(rhs)
    }

    fn zero() -> Self {
        Self::default()
    }

    fn one() -> Self {
        let mut coeffs = vec![F::ZERO; D];
        coeffs[0] = F::ONE;
        Self { coeffs }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c.is_zero())
    }

    fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c.neg()).collect();
        Self { coeffs }
    }

    fn scalar_mul(&self, scalar: F) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c * scalar).collect();
        Self { coeffs }
    }

    fn random() -> Self {
        Self::get_random()
    }
}

#[cfg(test)]
mod test {
    use crate::arith::{
        cyclotomic_ring::{CyclotomicRing, NTT_CUTOFF},
        field::Field64,
        ring::Ring,
    };
    use ark_ff::{AdditiveGroup, Field};
    use rand::{Rng, rng};

    #[test]
    /// Tests that the monomial multiplication is coherent with multiplication.
    fn test_monomial_mult() {
        for _ in 0..100 {
            const D: usize = 2;
            let mut monomial_coeffs = vec![Field64::ZERO; D];
            let monomial_non_null_term = rng().random_range(0..2 * D);

            if monomial_non_null_term < D {
                monomial_coeffs[monomial_non_null_term] = Field64::ONE;
            } else {
                monomial_coeffs[monomial_non_null_term % D] = -Field64::ONE;
            }

            let monomial = CyclotomicRing::<D, Field64> {
                coeffs: monomial_coeffs,
            };

            let polynomial = CyclotomicRing::get_random();

            let res_mul = polynomial.mul(&monomial);
            let res_monomial_mul = polynomial.multiply_by_monomial(monomial_non_null_term);

            assert_eq!(res_mul.coeffs, res_monomial_mul.coeffs);
        }
    }

    #[test]
    fn test_basic_mul() {
        const D: usize = 4;
        let lhs = CyclotomicRing::<D, Field64>::from_coeffs(&[
            Field64::from(1u64),
            Field64::from(2u64),
            Field64::from(0u64),
            Field64::from(3u64),
        ]);

        let rhs = CyclotomicRing::<D, Field64>::from_coeffs(&[
            Field64::from(0u64),
            Field64::from(1u64),
            Field64::from(9u64),
            Field64::from(0u64),
        ]);

        let res = lhs.basic_mul(&rhs);
        assert_eq!(
            res.coeffs,
            vec![
                Field64::from(18446744069414584318u64),
                Field64::from(18446744069414584295u64),
                Field64::from(11u64),
                Field64::from(18u64)
            ]
        );
    }

    #[test]
    fn test_ntt_mul() {
        const D: usize = NTT_CUTOFF * 2;
        let lhs = CyclotomicRing::<D, Field64>::get_random();
        let rhs = CyclotomicRing::<D, Field64>::get_random();

        let basic_res = lhs.basic_mul(&rhs);
        let ntt_res = lhs.try_ntt_mul(&rhs);
        assert_eq!(basic_res.coeffs, ntt_res.coeffs);
    }

    #[test]
    fn test_mul_via_ntt() {
        const D: usize = NTT_CUTOFF * 2;
        let lhs = CyclotomicRing::<D, Field64>::get_random();
        let rhs = CyclotomicRing::<D, Field64>::get_random();

        let basic_res = lhs.basic_mul(&rhs);
        let mul_res = lhs.mul(&rhs);
        assert_eq!(basic_res.coeffs, mul_res.coeffs);
    }

    #[test]
    fn test_mul_via_basic() {
        const D: usize = NTT_CUTOFF / 2;
        let lhs = CyclotomicRing::<D, Field64>::get_random();
        let rhs = CyclotomicRing::<D, Field64>::get_random();

        let basic_res = lhs.basic_mul(&rhs);
        let mul_res = lhs.mul(&rhs);
        assert_eq!(basic_res.coeffs, mul_res.coeffs);
    }
}
