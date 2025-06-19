use std::{iter::Sum, marker::PhantomData};

use ark_ff::PrimeField;
use ark_serialize::SerializationError;
use ark_std::rand::Rng;

use crate::{
    arith::ntt::{ntt::Ntt, tfhe_based_ntt::TfheBasedNtt},
    rand::get_rng,
};

use super::ring::Ring;

// At what point we start using NTT for multiplication.
const NTT_CUTOFF: usize = 32;

/// Represents an element of F\[X\]/(X^N + 1).
#[derive(Clone, Debug)]
pub struct CyclotomicRing<const D: usize, F: PrimeField, N: Ntt + Clone + Default = TfheBasedNtt> {
    pub coeffs: Vec<F>,
    _ntt: PhantomData<N>,
}

impl<const D: usize, F: PrimeField, N: Ntt + Clone + Default> CyclotomicRing<D, F, N> {
    pub fn from_coeffs(coeffs: &[F]) -> Self {
        assert!(
            coeffs.len() == D,
            "Coefficient vector must have length equal to the degree of the cyclotomic ring"
        );

        Self {
            coeffs: coeffs.to_vec(),
            _ntt: PhantomData,
        }
    }

    pub fn add_constant(&self, constant: F) -> Self {
        let mut res: Self = self.clone();
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

        Self {
            coeffs: rotated_coeffs,
            _ntt: PhantomData,
        }
    }

    pub fn get_random() -> Self {
        let mut rng = get_rng();

        let coeffs = (0..D).map(|_| F::rand(&mut rng)).collect();

        Self {
            coeffs,
            _ntt: PhantomData,
        }
    }

    pub fn get_random_bin() -> Self {
        let mut rng = get_rng();

        let coeffs = (0..D).map(|_| F::from(rng.gen_range(0..=1))).collect();

        Self {
            coeffs,
            _ntt: PhantomData,
        }
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
        Self {
            coeffs,
            _ntt: PhantomData,
        }
    }

    pub fn ntt_mul(&self, rhs: &Self) -> Self {
        Self {
            coeffs: N::mul::<D, F>(&self.coeffs, &rhs.coeffs),
            _ntt: PhantomData,
        }
    }

    pub fn serialize(&self, bytes: &mut Vec<u8>) -> Result<(), SerializationError> {
        for coeff in &self.coeffs {
            coeff.serialize_uncompressed(&mut *bytes)?;
        }
        Ok(())
    }
}

impl<const D: usize, F: PrimeField, N: Ntt + Clone + Default> Default for CyclotomicRing<D, F, N> {
    fn default() -> Self {
        Self {
            coeffs: vec![F::ZERO; D],
            _ntt: PhantomData,
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
        self.ntt_mul(rhs)
    }

    fn zero() -> Self {
        Self::default()
    }

    fn one() -> Self {
        let mut coeffs = vec![F::ZERO; D];
        coeffs[0] = F::ONE;
        Self {
            coeffs,
            _ntt: PhantomData,
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c.is_zero())
    }

    fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c.neg()).collect();
        Self {
            coeffs,
            _ntt: PhantomData,
        }
    }

    fn scalar_mul(&self, scalar: F) -> Self {
        let coeffs = self.coeffs.iter().map(|&c| c * scalar).collect();
        Self {
            coeffs,
            _ntt: PhantomData,
        }
    }

    fn random() -> Self {
        Self::get_random()
    }
}

#[cfg(test)]
mod test {
    use std::marker::PhantomData;

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
                _ntt: PhantomData,
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
        let ntt_res = lhs.ntt_mul(&rhs);
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
