fn main() {
    // println!("Hello world!")
}

#[cfg(test)]
mod test {
    use matvec::{
        arith::{
            cyclotomic_ring::CyclotomicRing,
            field::{Field64, Field64_2},
            linalg::Matrix,
        },
        protocol::{prover::Prover, verifier::Verifier},
        rlwe::{decrypt, encrypt},
    };

    #[test]
    fn test_functionality() {
        pub const D: usize = 1024;
        pub const INTEGER_WIDTH: usize = D * D;
        pub const INTEGER_HEIGHT: usize = 32;
        pub type F = Field64;
        pub type F2 = Field64_2;

        // generate matrix data
        let data = (1..=INTEGER_HEIGHT * INTEGER_WIDTH)
            .map(|x| F::from(x as u64))
            .collect::<Vec<_>>();

        let m = Matrix::from_vec(data, INTEGER_WIDTH);
        // println!("m: {:#?}", m);

        // generate vector data
        let v = (1..=INTEGER_WIDTH)
            .map(|x| (x % 16) as u64)
            .collect::<Vec<_>>();

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

        let (m_rq, m_polyring, m_mle, m_mle_over_base_f, z1_num_vars, mut transcript) =
            Prover::<D, F2>::preprocess(&m);

        let proof = Prover::<D, F2>::prove(
            &m_rq,
            &m_polyring,
            &m_mle,
            &m_mle_over_base_f,
            z1_num_vars,
            &mut transcript,
            &x,
        );

        let res = &proof
            .y
            .iter()
            .map(|c| decrypt(&sk, c)[0])
            .collect::<Vec<_>>();

        assert_eq!(m_v, *res);

        let (m_rq, z1_num_vars, mut transcript) = Verifier::<D, F2>::preprocess(&m);

        let verifier_res =
            Verifier::<D, F2>::verify(&m_rq, z1_num_vars, &mut transcript, &x, proof.clone());

        assert!(verifier_res.is_ok());
    }
}
