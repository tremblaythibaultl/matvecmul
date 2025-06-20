use std::iter::Sum;

use ark_ff::{Field, PrimeField};
use ark_std::rand::Rng;

use super::ring::Ring;

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

    // TODO: use NTT for better performances
    fn mul(&self, rhs: &CyclotomicRing<D, F>) -> Self {
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
    use crate::arith::{cyclotomic_ring::CyclotomicRing, field::Field64, ring::Ring};
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
}
