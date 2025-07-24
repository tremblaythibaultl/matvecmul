use ark_ff::{Field, PrimeField};

use crate::arith::cyclotomic_ring::CyclotomicRing;

use super::ring::Ring;

/// Represents an element of F\[X\] of maximal degree 2*D.
#[derive(Clone, Debug, Default)]
pub struct PolynomialRing<const D: usize, F: Field> {
    pub coeffs: Vec<F>,
}

impl<const D: usize, F: PrimeField> PolynomialRing<D, F> {
    pub fn long_division_by_cyclotomic(&mut self) -> (PolynomialRing<D, F>, CyclotomicRing<D, F>) {
        let mut quotient = PolynomialRing::<D, F>::zero();

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

                self.coeffs[quotient_term_degree] -= reduced_coeff;
                self.coeffs[i] = F::zero();
            }
        }

        // Get the cyclotomic ring element (i.e. `self` reduced modulo `X^D + 1`)
        let remainder = CyclotomicRing::<D, F>::from_coeffs(&self.coeffs[..D]);

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
        todo!()
    }

    fn sub(&self, other: &Self) -> Self {
        todo!()
    }

    fn sub_assign(&mut self, other: &Self) {
        todo!()
    }

    fn mul(&self, other: &Self) -> Self {
        todo!()
    }

    fn zero() -> Self {
        Self {
            coeffs: vec![F::zero(); D],
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
    use crate::arith::{field::Field64, polynomial_ring::PolynomialRing};

    #[test]
    fn test_long_division_by_cyclotomic() {
        const D: usize = 4;
        // let mut poly = PolynomialRing::<D, Field64>::random();
        let mut poly = PolynomialRing::<D, Field64> {
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
        assert_eq!(quotient.coeffs.len(), 4);
        assert_eq!(remainder.coeffs.len(), 4);

        assert_eq!(
            quotient.coeffs,
            vec![
                Field64::from(0),
                Field64::from(3),
                Field64::from(1),
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
}
