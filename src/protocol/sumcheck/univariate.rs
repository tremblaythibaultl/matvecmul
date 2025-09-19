use ark_ff::Field;

use crate::arith::linalg::gaussian_elimination;

/// Represents a univariate polynomial stored in coefficient form.
/// TODO: Check if coefficient form is the best way to encode the message polynomials
#[derive(Clone, Debug)]
pub struct UnivariatePolynomial<F> {
    coeffs: Vec<F>,
}

impl<F: Field> UnivariatePolynomial<F> {
    pub fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs }
    }

    pub fn coeffs(&self) -> &[F] {
        &self.coeffs
    }

    // TODO: double-check this is the best way to do polynomial interpolation.
    pub fn new_from_eval_points(evals: Vec<F>) -> Self {
        let n = evals.len();
        let xs = (0..n).map(|x| F::from(x as u64)).collect::<Vec<F>>();

        let mut vandermonde: Vec<Vec<F>> = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(n);
            let x = xs[i];
            row.push(F::one());
            row.push(x);
            for j in 2..n {
                row.push(row[j - 1] * x);
            }
            row.push(evals[i]);
            vandermonde.push(row);
        }

        Self {
            coeffs: gaussian_elimination(&mut vandermonde),
        }
    }

    pub fn eval(&self, evaluation_point: &F) -> F {
        let eval = (0..self.coeffs.len())
            .scan(F::ONE, |state, i| {
                let current = *state;
                *state *= evaluation_point;
                Some(current * self.coeffs[i])
            })
            .sum();

        eval
    }
}
