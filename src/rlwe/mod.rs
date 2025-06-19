use crate::arith::cyclotomic_ring::CyclotomicRing;
use ark_ff::PrimeField;

// MLWE rank
const K: usize = 1;
// Log of message modulus
const lg_p: usize = 4;

pub struct RLWE<const D: usize, F: PrimeField> {
    mask: Vec<CyclotomicRing<D, F>>,
    body: CyclotomicRing<D, F>,
}

// Performs encryption
pub fn encrypt<const D: usize, F: PrimeField>(
    sk: &Vec<CyclotomicRing<D, F>>,
    msg: &Vec<u64>,
) -> RLWE<D, F> {
    // ensure message is in message space
    for m in msg {
        assert!(m >> lg_p == 0);
    }

    // sample random mask
    let a: Vec<CyclotomicRing<D, F>> = (0..K).map(|_| CyclotomicRing::get_random()).collect();

    // compute body
    let b: CyclotomicRing<D, F> = a.iter().zip(sk).map(|(ai, si)| ai.mul(&si)).sum();

    // compute scaling factor (not exact)
    let delta = F::MODULUS_BIT_SIZE as u64 - lg_p as u64;

    let mut mu_coeffs = Vec::<F>::with_capacity(D);
    for i in 0..D {
        mu_coeffs.push(F::from(delta * msg[i]));
        // TODO: add error
    }

    RLWE {
        mask: a,
        body: b.add(&CyclotomicRing { coeffs: mu_coeffs }),
    }
}

pub fn decrypt<const D: usize, F: PrimeField>(
    sk: &Vec<CyclotomicRing<D, F>>,
    c: &RLWE<D, F>,
) -> Vec<F> {
    let b: CyclotomicRing<D, F> = c.mask.iter().zip(sk).map(|(ai, si)| ai.mul(&si)).sum();

    let mu_star = c.body.sub(&b);

    // TODO: round out error

    let mut res = Vec::<F>::with_capacity(D);

    let delta = F::MODULUS_BIT_SIZE as u64 - lg_p as u64;

    for i in 0..D {
        res.push(mu_star.coeffs[i] / F::from(delta)); // scale down value - not quite correct.
    }

    res
}

impl<const D: usize, F: PrimeField> RLWE<D, F> {
    pub fn mul_constant(&self, rhs: &CyclotomicRing<D, F>) -> RLWE<D, F> {
        let res_mask = self
            .mask
            .iter()
            .map(|p| p.mul(rhs))
            .collect::<Vec<CyclotomicRing<D, F>>>();
        let res_body = self.body.mul(rhs);

        RLWE {
            mask: res_mask,
            body: res_body,
        }
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
