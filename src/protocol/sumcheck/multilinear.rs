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

#[cfg(test)]
mod tests {

    use crate::{
        arith::field::Field64_2,
        protocol::{sumcheck::multilinear::MultilinearPolynomial, utils::build_eq_poly},
    };

    pub type F = Field64_2;

    #[test]
    fn test_fast_eval_structured() {
        let evals = vec![F::from(1), F::from(2), F::from(3), F::from(4)];
        let padding_factor = 8u32;
        let num_vars = 2;
        let padded_num_vars = num_vars + padding_factor.ilog2() as usize;

        let unpadded_mle = MultilinearPolynomial::new(evals, num_vars);

        let padded_mle_evals = (0..padding_factor)
            .map(|_| unpadded_mle.evals().iter().cloned())
            .flatten()
            .collect();

        let mut padded_mle = MultilinearPolynomial::new(padded_mle_evals, padded_num_vars);
        let challenges = vec![F::from(2), F::from(3), F::from(5), F::from(7), F::from(11)];

        let challenge_subvec = &challenges
            [(padding_factor.ilog2() as usize)..num_vars + (padding_factor.ilog2() as usize)];
        let mut eq_r_mle_evals = Vec::<F>::with_capacity(1 << num_vars);
        build_eq_poly(
            &challenge_subvec.iter().rev().cloned().collect::<Vec<F>>(),
            &mut eq_r_mle_evals,
        );

        let eval = unpadded_mle
            .evals()
            .iter()
            .zip(eq_r_mle_evals.iter())
            .map(|(a, b)| *a * *b)
            .sum::<F>();

        for chal in challenges {
            padded_mle.bind_to_challenge(&chal);
        }

        assert_eq!(eval, padded_mle.evals()[0]);
    }
}
