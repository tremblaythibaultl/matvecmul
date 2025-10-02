use ark_ff::Field;
use rayon::prelude::*;

use crate::protocol::{
    sample_random_challenge,
    sumcheck::{multilinear::MultilinearPolynomial, univariate::UnivariatePolynomial},
};

pub mod multilinear;
pub mod univariate;
#[derive(Clone)]
pub struct SumCheckProof<F: Field> {
    claim: F,
    polys: Vec<UnivariatePolynomial<F>>,
}

impl<F: Field> SumCheckProof<F> {
    pub fn verify(&self, num_vars: usize, max_degree: usize) -> Result<(F, F, Vec<F>), usize> {
        assert_eq!(self.polys.len(), num_vars);

        let mut claim = self.claim;
        let mut challenges = Vec::<F>::with_capacity(num_vars);

        for i in 0..self.polys.len() {
            // verify degree
            assert!(self.polys[i].coeffs().len() <= max_degree + 1);

            // TODO: implement fiat-shamir logic
            let challenge = sample_random_challenge::<F>(true);

            let zero_one_eval =
                self.polys[i].coeffs()[0] + self.polys[i].coeffs().iter().sum::<F>();

            let bound_poly = self.polys[i].eval(&challenge);

            if !(claim == zero_one_eval) {
                return Err(i);
            }

            claim = bound_poly;
            challenges.push(challenge);
        }

        return Ok((self.claim, claim, challenges));
    }
}

// Method to generate prover messages for sumcheck protocol for arbitrary-degree product of MLEs
// Inspired from Libra's sum-check for product of multilinear polynomials (ia.cr/2019/317).
pub fn prove<F: Field>(
    claim: F,
    mles: &mut Vec<MultilinearPolynomial<F>>,
    num_vars: usize,
) -> (SumCheckProof<F>, Vec<F>) {
    let max_degree = mles.len(); // Degree of the polynomial is equal to the number of MLEs being multiplied. We assume the factor of each MLE is `1`.

    let mut challenges = Vec::<F>::with_capacity(num_vars);
    let mut round_polys = Vec::<UnivariatePolynomial<F>>::with_capacity(num_vars);

    for _ in 0..num_vars {
        // we assume all MLEs have the same size
        let mle_half = mles[0].evals().len() / 2;

        let evals_to_sum: Vec<Vec<F>> = (0..mle_half)
            .into_par_iter()
            .map(|poly_term_i| {
                let mut evals = Vec::<F>::with_capacity(max_degree + 1);

                // term for t=0
                let t_zero_evals = mles
                    .iter()
                    .map(|poly| poly.evals()[poly_term_i])
                    .collect::<Vec<_>>();
                evals.push(t_zero_evals.iter().product());

                // term for t=1
                let t_one_evals = mles
                    .iter()
                    .map(|poly| poly.evals()[mle_half + poly_term_i])
                    .collect::<Vec<_>>();
                evals.push(t_one_evals.iter().product());

                let mut t_i_evals = t_one_evals;

                // terms for t=(2..=max_degree)
                for _ in 2..=max_degree {
                    let mut poly_evals = vec![F::zero(); mles.len()];
                    mles.iter().enumerate().for_each(|(i, mle)| {
                        poly_evals[i] = t_i_evals[i] + mle.evals()[mle_half + poly_term_i]
                            - mle.evals()[poly_term_i];
                    });
                    evals.push(poly_evals.iter().product::<F>());
                    t_i_evals = poly_evals;
                }
                evals
            })
            .collect();

        let round_poly_evaluations = evals_to_sum.into_par_iter().reduce(
            || vec![F::zero(); max_degree + 1],
            |mut acc, evals| {
                for (i, &eval) in evals.iter().enumerate() {
                    acc[i] += eval;
                }
                acc
            },
        );

        round_polys.push(UnivariatePolynomial::new_from_eval_points(
            round_poly_evaluations,
        ));

        let random_challenge = sample_random_challenge::<F>(true);
        challenges.push(random_challenge); // TODO: Implement Fiat-Shamir logic

        // bind to verifier's challenge
        mles.iter_mut()
            .for_each(|mle| mle.bind_to_challenge(&random_challenge));
    }

    (
        SumCheckProof {
            claim,
            polys: round_polys,
        },
        challenges,
    )
}

#[cfg(test)]
mod tests {
    use crate::{
        arith::field::Field64_2,
        protocol::sumcheck::{multilinear::MultilinearPolynomial, prove},
    };

    pub type F = Field64_2;

    #[test]
    fn test_sum_check() {
        let evals_1 = vec![
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4),
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4),
        ];
        let evals_2 = vec![
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
        ];

        let claim = evals_1
            .iter()
            .zip(evals_2.iter())
            .map(|(a, b)| a * b)
            .sum::<F>();

        println!("Claim: {:?}", claim);

        let num_vars = 3;
        let max_degree = 2;

        let mle_1 = MultilinearPolynomial::new(evals_1, num_vars);
        let mle_2 = MultilinearPolynomial::new(evals_2, num_vars);

        let (proof, _challenges) = prove(claim, &mut vec![mle_1, mle_2], num_vars);

        let v = proof.verify(num_vars, max_degree);

        println!("Verification result: {:?}", v);
    }
}
