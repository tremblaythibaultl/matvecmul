use core::num;

use ark_ff::Field;

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, linalg::Matrix, polynomial_ring::PolynomialRing},
    protocol::{
        sumcheck::{SumCheckProof, multilinear::MultilinearPolynomial},
        utils::build_eq_poly,
    },
    rlwe::RLWE,
};

pub mod prover;
pub mod sumcheck;
mod utils;
pub mod verifier;
#[derive(Clone)]
pub struct Proof<const D: usize, F: Field> {
    pub y: Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    pub r: Vec<Vec<PolynomialRing<D, F::BasePrimeField>>>, // should be a commitment to `r`
    pub z1_sumcheck_proof: SumCheckProof<F>,
}

pub fn sample_random_challenge<F: Field>(test_mode: bool) -> F {
    if test_mode {
        // In test mode, return a fixed challenge

        // let coeffs = (1..=F::extension_degree())
        //     .map(|i| F::BasePrimeField::from(i as u64))
        //     .collect::<Vec<_>>();

        // return F::from_base_prime_field_elems(coeffs).unwrap(); // Should be safe to unwrap as the slice length matches the extension degree
        return F::from(11);
    } else {
        // TODO: implement a proper random challenge sampling
        // Should basically hash everything that is known to the verifier, i.e. `M`, `x`, `y`, and `r`.
        // The random challenge will probably be an element of `Gal(p, d)` for some `d` (probably 2?) so we should choose an appropriate hash function (Poseidon?)
        todo!()
    }
}

pub fn preprocess<const D: usize, F: Field>(
    m: &Matrix<F::BasePrimeField>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
) -> (
    Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    Matrix<PolynomialRing<D, F::BasePrimeField>>,
    Vec<RLWE<PolynomialRing<D, F::BasePrimeField>>>,
    MultilinearPolynomial<F>,
    usize,
) {
    // interpret each row as a ring elements and rearrange columns
    let m_rq = m.lift_to_rq::<D>();

    // lift matrix entries to polynomial ring elements of max degree 2*D
    // this performs a clone of the coefficients
    // should look into optimizing
    let m_polyring = m_rq.lift_to_polynomial_ring();

    let x_polyring = x
        .iter()
        .map(|elem| elem.lift_to_polynomial_ring())
        .collect::<Vec<_>>();

    // this also clones the coefficients. should look into optimizing.
    let (m_mle_evals, num_vars): (Vec<F::BasePrimeField>, usize) = m_rq.to_mle_evals();

    let m_mle = MultilinearPolynomial::new(
        m_mle_evals
            .into_iter()
            .map(|e| F::from_base_prime_field(e))
            .collect(),
        num_vars,
    );

    (m_rq, m_polyring, x_polyring, m_mle, num_vars)
}

pub fn compute_mles<const D: usize, F: Field>(
    m_rq: &Matrix<CyclotomicRing<D, F::BasePrimeField>>,
    m_mle: &MultilinearPolynomial<F>,
    x: &Vec<RLWE<CyclotomicRing<D, F::BasePrimeField>>>,
    num_vars: usize,
    alpha: &F,
    tau: &F,
) -> Vec<MultilinearPolynomial<F>> {
    let powers_of_alpha: Vec<F> = (0..D)
        .scan(F::one(), |state, _| {
            let result = *state;
            *state *= alpha;
            Some(result)
        })
        .collect();

    // this might not be very efficient memory-wise, but we need it to keep genericity in the sumcheck prover.
    // TODO: send minimal information to the sumcheck prover and a function describing how to pad the MLEs
    let padded_alpha_mle_evals = (0..m_rq.width() * m_rq.height())
        .flat_map(|_| powers_of_alpha.clone())
        .collect::<Vec<_>>();

    let alpha_mle = MultilinearPolynomial::new(padded_alpha_mle_evals, num_vars);

    // only consider the mask for now. will probably need to run the protocol K+1 times where K is the RLWE rank
    let x_alpha_mle_evals = x
        .iter()
        .map(|ct| {
            ct.get_ring_elements()[0]
                .coeffs
                .iter()
                .zip(powers_of_alpha.iter())
                .map(|(c, a)| a.mul_by_base_prime_field(c))
                .sum::<F>()
        })
        .collect::<Vec<F>>();

    // pad x_alpha_mle_evals
    let padded_x_alpha_mle_evals = (0..m_rq.height() as usize)
        .flat_map(|_| {
            x_alpha_mle_evals
                .iter()
                .map(|e| vec![*e; D])
                .flatten()
                .collect::<Vec<F>>()
        })
        .collect::<Vec<F>>();

    let x_alpha_mle: MultilinearPolynomial<F> =
        MultilinearPolynomial::new(padded_x_alpha_mle_evals, num_vars);

    // compute eq polynomial
    let mut eq_tau_mle_evals = Vec::<F>::with_capacity(num_vars);
    let vec_tau = vec![*tau; num_vars];
    build_eq_poly(&vec_tau, &mut eq_tau_mle_evals);

    let eq_tau_mle = MultilinearPolynomial::new(eq_tau_mle_evals, num_vars);

    vec![
        eq_tau_mle,
        m_mle.clone(),
        x_alpha_mle.clone(),
        alpha_mle.clone(),
    ]
}
