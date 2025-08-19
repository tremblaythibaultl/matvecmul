use ark_ff::Field;

/// Represents a multilinear polynomial stored in evaluation form over {0, 1}^num_variables according to lexicographic ordering.
pub struct MultilinearPolynomial<F> {
    num_variables: usize,
    evals: Vec<F>,
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn new(evals: Vec<F>, num_variables: usize) -> Self {
        Self {
            num_variables,
            evals,
        }
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn sum_over_boolean_hypercube(&self) -> F {
        (0..(1 << self.num_variables))
            .map(|point| self.evals[point])
            .sum()
    }
}
