use ark_ff::{Field, PrimeField};

use crate::arith::cyclotomic_ring::CyclotomicRing;

use super::ring::Ring;

/// Represents an element of F\[X\] of maximal degree D.
// We keep a constant degree D for the polynomial ring because the elements of interest are polynomials of degree at most `2*D_cyclo`, where `D_cyclo` is the degree of the cyclotomic ring.
#[derive(Clone, Debug, Default)]
pub struct PolynomialRing<const D: usize, F: Field> {
    pub coeffs: Vec<F>,
}

impl<const D: usize, F: PrimeField> PolynomialRing<D, F> {
    pub fn long_division_by_cyclotomic<const D_cyclo: usize>(
        &mut self,
    ) -> (PolynomialRing<D_cyclo, F>, CyclotomicRing<D_cyclo, F>) {
        assert!(
            D_cyclo * 2 == D,
            "the cyclotomic degree must be half the maximum polynomial degree"
        );

        let mut quotient = PolynomialRing::<D_cyclo, F>::zero();

        // Perform long division of `self` by the cyclotomic polynomial `X^{D_cyclo} + 1`
        // Only need to iterate through the last `D_cyclo` coefficients of `self` because `self` is of degree at most `2 * D_cyclo`.
        for i in (D_cyclo..2 * D_cyclo).rev() {
            if self.coeffs[i] != F::zero() {
                let quotient_term_degree = i - D_cyclo;

                // The coefficient of `self` that corresponds to `X^{i}`
                let reduced_coeff = self.coeffs[i];

                // Add `reduced_coeff * X^{quotient_term_degree}` to the quotient
                // `quotient_term_degree` is monotonically decreasing, so we can safely assign the coefficient
                quotient.coeffs[quotient_term_degree] = reduced_coeff;

                self.coeffs[quotient_term_degree] -= reduced_coeff;
                self.coeffs[i] = F::zero();
            }
        }

        // Get the cyclotomic ring element (i.e. `self` reduced modulo `X^{D_cyclo} + 1`)
        let remainder = CyclotomicRing::<D_cyclo, F>::from_coeffs(&self.coeffs[..D_cyclo]);

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
        const D: usize = 8;
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
        let (quotient, remainder) = poly.long_division_by_cyclotomic::<4>();
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
