use ark_ff::Field;

use crate::protocol::sumcheck::multilinear::MultilinearPolynomial;

pub fn build_eq_poly<F: Field>(r: &[F], buf: &mut Vec<F>) {
    if r.len() == 1 {
        buf.push(F::one() - r[0]);
        buf.push(r[0]);
    } else {
        build_eq_poly(&r[1..], buf);

        let mut res = vec![F::zero(); buf.len() << 1];
        res.iter_mut().enumerate().for_each(|(i, val)| {
            let bi = buf[i >> 1];
            let tmp = r[0] * bi;
            if i & 1 == 0 {
                *val = bi - tmp;
            } else {
                *val = tmp;
            }
        });
        *buf = res;
    }
}

// Computes the sum over the boolean hypercube of the product of a vector of MLEs
pub fn sum_over_boolean_hypercube<F: Field>(mles: &[MultilinearPolynomial<F>]) -> F {
    let len = mles[0].evals().len();

    // Check all MLEs have the same length
    for mle in mles.iter() {
        assert_eq!(mle.evals().len(), len);
    }

    let mut sum = F::ZERO;

    for i in 0..len {
        let prod: F = mles.iter().map(|mle| mle.evals()[i]).product();
        sum += prod;
    }

    sum
}
