use ark_ff::PrimeField;

pub trait Ntt {
    fn mul<const D: usize, F: PrimeField>(lhs: &[F], rhs: &[F]) -> Vec<F>;
}
