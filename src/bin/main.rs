fn main() {}

#[cfg(test)]
mod test {
    use ark_ff::PrimeField;
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
        pub const D: usize = 1 << 10;
        pub const P: usize = 1 << 4;
        pub const INTEGER_WIDTH: usize = 1 << 20;
        pub const INTEGER_HEIGHT: usize = 1 << 6;
        pub type F = Field64;
        pub type F2 = Field64_2;

        // generate matrix data
        let data = (1..=INTEGER_HEIGHT * INTEGER_WIDTH)
            .map(|x| F::from((x % P) as u64))
            .collect::<Vec<_>>();

        let m = Matrix::from_vec(data, INTEGER_WIDTH);

        // generate vector data
        let v = (1..=INTEGER_WIDTH)
            .map(|x| (x % P) as u64)
            .collect::<Vec<_>>();

        // compute plaintext matrix-vector multiplication
        let m_v = m
            .mat_vec_mul(&v.iter().map(|&v_i| F::from(v_i)).collect::<Vec<F>>())
            .iter()
            .map(|y_i| {
                let y_i_u64 = y_i.into_bigint().0[0];
                F::from(y_i_u64 % P as u64)
            })
            .collect::<Vec<F>>();

        // generate a secret key
        let sk = vec![CyclotomicRing::get_random_bin()];

        // encrypt the vector
        let x = v
            .chunks(D)
            .map(|subvec| encrypt::<D, P, F>(&sk, &subvec))
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
            true,
        );

        let res = &proof
            .y
            .iter()
            .map(|c| decrypt::<D, P, F>(&sk, c)[0])
            .collect::<Vec<_>>();

        assert_eq!(m_v, *res);

        let (m_rq, z1_num_vars, mut transcript) = Verifier::<D, F2>::preprocess(&m);

        let verifier_res =
            Verifier::<D, F2>::verify(&m_rq, z1_num_vars, &mut transcript, &x, proof.clone());

        assert!(verifier_res.is_ok());
    }
}
