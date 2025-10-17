use ark_ff::Field;

use crate::rand::get_rng;

pub trait Ring: Clone + Default {
    /// The degree of the ring extension over the base field
    const DEGREE: usize;

    /// The base field type
    type BaseField: Field;

    fn add(&self, other: &Self) -> Self;

    fn add_assign(&mut self, other: &Self);

    fn sub(&self, other: &Self) -> Self;

    fn sub_assign(&mut self, other: &Self);

    fn mul(&self, other: &Self) -> Self;

    fn zero() -> Self;

    fn one() -> Self;

    fn is_zero(&self) -> bool;

    fn neg(&self) -> Self;

    /// Scalar multiplication by base field element
    fn scalar_mul(&self, scalar: Self::BaseField) -> Self;

    /// Get random element
    fn random() -> Self;
}

impl<F: Field> Ring for F {
    const DEGREE: usize = 1;

    type BaseField = F;

    fn add(&self, other: &Self) -> Self {
        *self + *other
    }

    fn add_assign(&mut self, other: &Self) {
        *self += *other;
    }

    fn sub(&self, other: &Self) -> Self {
        *self - *other
    }

    fn sub_assign(&mut self, other: &Self) {
        *self -= *other;
    }

    fn mul(&self, other: &Self) -> Self {
        *self * *other
    }

    fn zero() -> Self {
        F::zero()
    }

    fn one() -> Self {
        F::one()
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn neg(&self) -> Self {
        -(*self)
    }

    fn scalar_mul(&self, scalar: Self::BaseField) -> Self {
        *self * scalar
    }

    fn random() -> Self {
        F::rand(&mut get_rng())
    }
}
