use ark_ff::{Field, PrimeField};

use crate::arith::cyclotomic_ring::CyclotomicRing;

use super::ring::Ring;

/// Represents an element of F\[X\] of maximal degree 2*D.
#[derive(Clone, Debug, Default)]
pub struct PolynomialRing<const D: usize, F: Field> {
    pub coeffs: Vec<F>,
}

impl<const D: usize, F: PrimeField> PolynomialRing<D, F> {
    pub fn from_cyclotomic(cyclotomic: &CyclotomicRing<D, F>) -> Self {
        assert_eq!(
            cyclotomic.coeffs.len(),
            D,
            "Cyclotomic polynomial must have degree D"
        );

        // Careful that this clone is not too expensive.
        let mut coeffs = cyclotomic.coeffs.clone();

        // pad to degree 2*D
        coeffs.resize(2 * D, F::zero());

        Self { coeffs }
    }

    pub fn long_division_by_cyclotomic(&self) -> (PolynomialRing<D, F>, CyclotomicRing<D, F>) {
        let mut quotient = PolynomialRing::<D, F>::zero();
        let mut coeffs = self.coeffs.clone();
        // Perform long division of `self` by the cyclotomic polynomial `X^{D} + 1`
        // Only need to iterate through the last `D` coefficients of `self` because `self` is of degree at most `2 * D`.
        for i in (D..2 * D).rev() {
            if self.coeffs[i] != F::zero() {
                let quotient_term_degree = i - D;

                // The coefficient of `self` that corresponds to `X^{i}`
                let reduced_coeff = self.coeffs[i];

                // Add `reduced_coeff * X^{quotient_term_degree}` to the quotient
                // `quotient_term_degree` is monotonically decreasing, so we can safely assign the coefficient
                quotient.coeffs[quotient_term_degree] = reduced_coeff;

                coeffs[quotient_term_degree] -= reduced_coeff;
                coeffs[i] = F::zero();
            }
        }

        // Get the cyclotomic ring element (i.e. `self` reduced modulo `X^D + 1`)
        let remainder = CyclotomicRing::<D, F>::from_coeffs(&coeffs[..D]);

        (quotient, remainder)
    }
}

impl<const D: usize, F: Field> Ring for PolynomialRing<D, F> {
    const DEGREE: usize = D;

    type BaseField = F;

    fn add(&self, other: &Self) -> Self {
        todo!()
    }

    fn add_assign(&mut self, other: &Self) {
        self.coeffs
            .iter_mut()
            .zip(other.coeffs.iter())
            .for_each(|(a, b)| *a += *b);
    }

    fn sub(&self, other: &Self) -> Self {
        todo!()
    }

    fn sub_assign(&mut self, other: &Self) {
        todo!()
    }

    fn mul(&self, other: &Self) -> Self {
        for i in D..2 * D {
            assert!(
                self.coeffs[i] == F::zero() && other.coeffs[i] == F::zero(),
                "PolynomialRing multiplicands must have degree at most D"
            );
        }

        let mut coeffs = vec![F::zero(); 2 * D];

        for i in 0..D {
            for j in 0..D {
                if i + j < 2 * D {
                    coeffs[i + j] += self.coeffs[i] * other.coeffs[j];
                }
            }
        }

        Self { coeffs }
    }

    fn zero() -> Self {
        Self {
            coeffs: vec![F::zero(); 2 * D],
        }
    }

    fn one() -> Self {
        todo!()
    }

    fn is_zero(&self) -> bool {
        todo!()
    }

    fn neg(&self) -> Self {
        todo!()
    }

    fn scalar_mul(&self, scalar: Self::BaseField) -> Self {
        todo!()
    }

    fn random() -> Self {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use crate::arith::{
        cyclotomic_ring::CyclotomicRing, field::Field64, polynomial_ring::PolynomialRing,
        ring::Ring,
    };
    const D: usize = 4;

    #[test]
    fn test_long_division_by_cyclotomic() {
        let poly = PolynomialRing::<D, Field64> {
            coeffs: vec![
                Field64::from(1),
                Field64::from(0),
                Field64::from(2),
                Field64::from(0),
                Field64::from(0),
                Field64::from(3),
                Field64::from(1),
                Field64::from(0),
            ],
        };
        let (quotient, remainder) = poly.long_division_by_cyclotomic();

        assert_eq!(
            quotient.coeffs,
            vec![
                Field64::from(0),
                Field64::from(3),
                Field64::from(1),
                Field64::from(0),
                Field64::from(0),
                Field64::from(0),
                Field64::from(0),
                Field64::from(0)
            ]
        );
        assert_eq!(
            remainder.coeffs,
            vec![
                Field64::from(1),
                Field64::from(-3),
                Field64::from(1),
                Field64::from(0)
            ]
        );
    }

    #[test]
    fn test_mul_and_long_division() {
        let mut coeffs1 = vec![
            Field64::from(1),
            Field64::from(2),
            Field64::from(3),
            Field64::from(4),
        ];

        let mut coeffs2 = vec![
            Field64::from(-5),
            Field64::from(-6),
            Field64::from(-7),
            Field64::from(-8),
        ];

        let cyclo1 = CyclotomicRing::<D, Field64>::from_coeffs(&coeffs1);

        coeffs1.resize(2 * D, Field64::from(0));

        let poly1 = PolynomialRing::<D, Field64> { coeffs: coeffs1 };

        let cyclo2 = CyclotomicRing::<D, Field64>::from_coeffs(&coeffs2);

        coeffs2.resize(2 * D, Field64::from(0));

        let poly2 = PolynomialRing::<D, Field64> { coeffs: coeffs2 };

        let poly3 = poly1.mul(&poly2);

        let cyclo3 = cyclo1.mul(&cyclo2);

        let (_, remainder) = poly3.long_division_by_cyclotomic();

        assert_eq!(remainder.coeffs, cyclo3.coeffs);
    }
}
