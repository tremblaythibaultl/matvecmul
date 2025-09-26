use ark_ff::Field;

/// Represents a multilinear polynomial stored in evaluation form over {0, 1}^num_variables according to lexicographic ordering.
#[derive(Clone, Debug)]
pub struct MultilinearPolynomial<F> {
    num_variables: usize,
    evals: Vec<F>,
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn new(evals: Vec<F>, num_variables: usize) -> Self {
        assert_eq!(evals.len(), 1 << num_variables);

        Self {
            num_variables,
            evals,
        }
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn num_variables(&self) -> usize {
        self.num_variables
    }

    pub fn evals_mut(&mut self) -> &mut [F] {
        &mut self.evals
    }

    pub fn sum_over_boolean_hypercube(&self) -> F {
        self.evals.iter().sum()
    }

    pub fn bind_to_challenge(&mut self, challenge: &F) {
        let mle_half = self.evals().len() / 2;
        let (left, right) = self.evals_mut().split_at_mut(mle_half);

        left.iter_mut().zip(right.iter()).for_each(|(a, b)| {
            *a += *challenge * (*b - *a);
        });

        self.num_variables -= 1;
        self.evals.truncate(mle_half);
    }
}
