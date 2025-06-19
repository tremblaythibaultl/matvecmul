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
        pub const WIDTH: usize = D;
        pub type F = Field64;

        // Use a simple, predictable matrix for testing
        let data = vec![
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4),
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8),
            F::from(9),
            F::from(10),
            F::from(11),
            F::from(12),
            F::from(13),
            F::from(14),
            F::from(15),
            F::from(16),
        ];

        let m = Matrix::from_vec(data, WIDTH);

        let v = vec![1, 2, 3, 4];

        let sk = vec![CyclotomicRing::get_random_bin()];

        let x = encrypt(&sk, &v);

        let proof = Prover::<D, F>::prove(m, x);

        let m0 = decrypt(&sk, &proof.y[0]);
        let m1 = decrypt(&sk, &proof.y[1]);
        let m2 = decrypt(&sk, &proof.y[2]);
        let m3 = decrypt(&sk, &proof.y[3]);

        println!("m0: {:?}", m0);
        println!("m1: {:?}", m1);
        println!("m2: {:?}", m2);
        println!("m3: {:?}", m3);
    }
}
