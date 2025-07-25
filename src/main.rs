mod arith;
mod protocol;
mod rlwe;

fn main() {
    // println!("Hello world!")
}

#[cfg(test)]
mod test {
    use crate::{
        arith::{cyclotomic_ring::CyclotomicRing, field::Field64, linalg::Matrix},
        protocol::prover::Prover,
        rlwe::{decrypt, encrypt},
    };

    #[test]
    fn test_functionality() {
        pub const D: usize = 4;
        pub const WIDTH: usize = 2 * D;
        pub const HEIGHT: usize = D;
        pub type F = Field64;

        // generate matrix data
        let data = (1..=HEIGHT * WIDTH)
            .map(|x| F::from(x as u64))
            .collect::<Vec<_>>();

        let m = Matrix::from_vec(data, WIDTH);
        // println!("m: {:#?}", m);

        // generate vector data
        let v = (1..=WIDTH).map(|x| x as u64).collect::<Vec<_>>();

        // compute plaintext matrix-vector multiplication
        let m_v = m.mat_vec_mul(&v.iter().map(|&x| F::from(x)).collect::<Vec<F>>());

        // generate a secret key
        let sk = vec![CyclotomicRing::get_random_bin()];

        // encrypt the vector
        let x = v
            .chunks(D)
            .map(|subvec| encrypt(&sk, &subvec))
            .collect::<Vec<_>>();

        // ask the prover to compute the encrypted matrix-vector multiplication and return the result
        let proof = Prover::<D, F>::prove(&m, &x);

        let res = &proof
            .y
            .iter()
            .map(|c| decrypt(&sk, c)[0])
            .collect::<Vec<_>>();

        assert_eq!(m_v, *res);
    }
}
