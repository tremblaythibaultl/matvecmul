use ark_ff::PrimeField;

pub trait Ntt {
    fn mul<F: PrimeField>(lhs: &[F], rhs: &[F]) -> Vec<F>;
}
