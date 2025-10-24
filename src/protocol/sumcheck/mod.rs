use ark_ff::Field;
use rayon::prelude::*;

use crate::protocol::{
    sumcheck::{multilinear::MultilinearPolynomial, univariate::UnivariatePolynomial},
    transcript::{self, Blake3Transcript},
};

pub mod multilinear;
pub mod univariate;
#[derive(Clone)]
pub struct SumCheckProof<F: Field> {
    claim: F,
    polys: Vec<UnivariatePolynomial<F>>,
}

impl<F: Field> SumCheckProof<F> {
    pub fn verify(
        &self,
        num_vars: usize,
        max_degree: usize,
        transcript: &mut Blake3Transcript<F>,
    ) -> Result<(F, F, Vec<F>), usize> {
        assert_eq!(self.polys.len(), num_vars);

        let mut claim = self.claim;
        let mut challenges = Vec::<F>::with_capacity(num_vars);

        transcript.absorb(&claim);

        for i in 0..self.polys.len() {
            // verify degree
            assert!(self.polys[i].coeffs().len() <= max_degree + 1);

            let challenge = transcript.squeeze(1)[0];

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
    transcript: &mut transcript::Blake3Transcript<F>,
) -> (SumCheckProof<F>, Vec<F>) {
    let max_degree = mles.len(); // Degree of the polynomial is equal to the number of MLEs being multiplied. We assume the factor of each MLE is `1`.

    let mut challenges = Vec::<F>::with_capacity(num_vars);
    let mut round_polys = Vec::<UnivariatePolynomial<F>>::with_capacity(num_vars);

    transcript.absorb(&claim);

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

        let random_challenge = transcript.squeeze(1)[0];

        challenges.push(random_challenge);

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
        protocol::{
            sumcheck::{multilinear::MultilinearPolynomial, prove},
            transcript,
        },
    };

    pub type F = Field64_2;

    #[test]
    fn test_sum_check() {
        let evals_1 = vec![
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4),
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4),
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
        ];
        let evals_2 = vec![
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
            F::from(9),
            F::from(10),
            F::from(11),
            F::from(12),
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
            F::from(9),
            F::from(10),
            F::from(11),
            F::from(12),
        ];

        let claim = evals_1
            .iter()
            .zip(evals_2.iter())
            .map(|(a, b)| a * b)
            .sum::<F>();

        let num_vars = 4;
        let max_degree = 2;

        let mut mle_1 = MultilinearPolynomial::new(evals_1, num_vars);
        let mut mle_2 = MultilinearPolynomial::new(evals_2, num_vars);

        let mut p_transcript = transcript::Blake3Transcript::<F>::new();

        let (proof, _) = prove(
            claim,
            &mut vec![mle_1.clone(), mle_2.clone()],
            num_vars,
            &mut p_transcript,
        );

        let mut v_transcript = transcript::Blake3Transcript::<F>::new();

        let (_, final_claim, verifier_challenges) = proof
            .verify(num_vars, max_degree, &mut v_transcript)
            .unwrap();

        for chal in verifier_challenges {
            mle_1.bind_to_challenge(&chal);
            mle_2.bind_to_challenge(&chal);
        }

        let mle_1_eval = mle_1.evals()[0];
        let mle_2_eval = mle_2.evals()[0];
        let final_eval = mle_1_eval * mle_2_eval;

        assert_eq!(final_eval, final_claim)
    }
}
