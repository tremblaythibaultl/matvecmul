use crate::arith::cyclotomic_ring::CyclotomicRing;
use crate::arith::polynomial_ring::PolynomialRing;
use crate::arith::ring::Ring;
use ark_ff::PrimeField;

// MLWE rank
pub const K: usize = 1;
// Log of message modulus
const LG_P: usize = 4;

#[derive(Debug, Clone)]
pub struct RLWE<R> {
    mask: Vec<R>,
    body: R,
}

// Performs encryption
pub fn encrypt<const D: usize, F: PrimeField>(
    sk: &Vec<CyclotomicRing<D, F>>,
    msg: &[u64],
) -> RLWE<CyclotomicRing<D, F>> {
    // ensure message is in message space
    for m in msg {
        assert!(
            m >> LG_P == 0,
            "Message value {} is too large for the given plaintext modulus {}",
            m,
            1 << LG_P
        );
    }

    // ensure the message is of the correct length
    assert!(
        msg.len() == D,
        "Message length {} does not match the cyclotomic ring degree {}",
        msg.len(),
        D
    );

    // sample random mask
    let a: Vec<CyclotomicRing<D, F>> = (0..K).map(|_| CyclotomicRing::get_random()).collect();

    // compute body
    let b: CyclotomicRing<D, F> = a.iter().zip(sk).map(|(ai, si)| ai.mul(&si)).sum();

    // compute scaling factor (not exact)
    let delta = F::MODULUS_BIT_SIZE as u64 - LG_P as u64;

    let mut mu_coeffs = Vec::<F>::with_capacity(D);
    for i in 0..D {
        mu_coeffs.push(F::from(delta * msg[i]));
        // TODO: add error
    }

    RLWE {
        mask: a,
        body: b.add(&CyclotomicRing::from_coeffs(&mu_coeffs)),
    }
}

pub fn decrypt<const D: usize, F: PrimeField>(
    sk: &Vec<CyclotomicRing<D, F>>,
    c: &RLWE<CyclotomicRing<D, F>>,
) -> Vec<F> {
    let b: CyclotomicRing<D, F> = c.mask.iter().zip(sk).map(|(ai, si)| ai.mul(&si)).sum();

    let mu_star = c.body.sub(&b);

    // TODO: round out error

    let mut res = Vec::<F>::with_capacity(D);

    let delta = F::MODULUS_BIT_SIZE as u64 - LG_P as u64;

    for i in 0..D {
        res.push(mu_star.coeffs[i] / F::from(delta)); // scale down value - not quite correct.
    }

    res
}

impl<const D: usize, F: PrimeField> RLWE<CyclotomicRing<D, F>> {
    pub fn lift_to_polynomial_ring(&self) -> RLWE<PolynomialRing<D, F>> {
        let mask = self
            .mask
            .iter()
            .map(|elem| PolynomialRing::<D, F>::from_cyclotomic(elem))
            .collect::<Vec<_>>()
            .into_iter()
            .collect();

        let body = PolynomialRing::<D, F>::from_cyclotomic(&self.body);

        RLWE { mask, body }
    }
}

impl<const D: usize, F: PrimeField> RLWE<PolynomialRing<D, F>> {
    pub fn long_division_by_cyclotomic(
        &self,
    ) -> (Vec<PolynomialRing<D, F>>, Vec<CyclotomicRing<D, F>>) {
        let cap = self.mask.len() + 1;
        let mut quotients = Vec::with_capacity(cap);
        let mut remainders = Vec::with_capacity(cap);

        for elem in &self.mask {
            let (quotient, remainder) = elem.long_division_by_cyclotomic();
            quotients.push(quotient);
            remainders.push(remainder);
        }

        let (quotient_body, remainder_body) = self.body.long_division_by_cyclotomic();

        quotients.push(quotient_body);
        remainders.push(remainder_body);

        (quotients, remainders)
    }
}

impl<R: Ring> RLWE<R> {
    pub fn zero() -> Self {
        RLWE {
            mask: vec![R::zero(); K],
            body: R::zero(),
        }
    }

    pub fn add_assign(&mut self, rhs: &RLWE<R>) {
        for (p1, p2) in self.mask.iter_mut().zip(&rhs.mask) {
            p1.add_assign(p2);
        }
        self.body.add_assign(&rhs.body);
    }

    pub fn mul_constant(&self, rhs: &R) -> RLWE<R> {
        let res_mask = self.mask.iter().map(|p| p.mul(rhs)).collect::<Vec<R>>();
        let res_body = self.body.mul(rhs);

        RLWE {
            mask: res_mask,
            body: res_body,
        }
    }

    // pub fn get_ring_elements_ref(self) -> &[R] {
    //     &self.mask
    // }

    pub fn get_ring_elements(&self) -> Vec<R> {
        let mut elements = self.mask.clone();
        elements.push(self.body.clone());
        elements
    }
}

#[cfg(test)]
mod test {
    use crate::arith::{cyclotomic_ring::CyclotomicRing, field::Field64};

    use super::{decrypt, encrypt};

    #[test]
    fn test_enc_dec() {
        const D: usize = 2;

        let msg_enc = vec![2, 3];

        let sk = vec![CyclotomicRing::<D, Field64>::get_random_bin()];

        let c = encrypt(&sk, &msg_enc);

        let msg_dec = decrypt(&sk, &c);

        assert_eq!(
            msg_dec,
            vec![Field64::from(msg_enc[0]), Field64::from(msg_enc[1])]
        )
    }
}
