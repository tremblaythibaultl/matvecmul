use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use matvec::{
    arith::{
        cyclotomic_ring::CyclotomicRing,
        field::{Field64, Field64_2},
        linalg::Matrix,
    },
    protocol::{
        pcs::whir::Whir, prover::Prover, sumcheck::multilinear::MultilinearPolynomial,
        verifier::Verifier,
    },
    rand::get_rng,
    rlwe::{RLWE, decrypt, encrypt},
};
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    slice::ParallelSlice,
};
use whir::crypto::fields::FieldWithSize;

pub const D: usize = 1 << 10;
pub const TS: &[usize] = &[
    (1 << 1),
    (1 << 2),
    (1 << 3),
    (1 << 4),
    (1 << 5),
    (1 << 6),
    (1 << 7),
    (1 << 8),
    (1 << 9),
    (1 << 10),
];
pub const HEIGHTS: &[usize] = &[(1 << 1), (1 << 2), (1 << 3), (1 << 4), (1 << 5), (1 << 6)];
pub type F = Field64;
pub type F2 = Field64_2;

fn setup_benchmark_data(
    t: usize,
    height: usize,
) -> (
    Matrix<F>,
    Vec<F>,
    Vec<u64>,
    Vec<CyclotomicRing<D, F>>,
    Vec<RLWE<CyclotomicRing<D, F>>>,
) {
    // Generate matrix data
    let integer_width = D * t;
    let integer_height = height;

    let data = (1..=integer_height * integer_width)
        .map(|x| F::from((x % 16) as u64))
        .collect::<Vec<_>>();
    let m = Matrix::from_vec(data, integer_width);

    // Generate vector data
    let v = (1..=integer_width)
        .map(|x| (x % 16) as u64)
        .collect::<Vec<_>>();

    // Compute plaintext matrix-vector multiplication for verification
    let m_v = m.mat_vec_mul(&v.iter().map(|&x| F::from(x)).collect::<Vec<F>>());

    // Generate a secret key
    let sk = vec![CyclotomicRing::get_random_bin()];

    let x = v
        .chunks(D)
        .map(|subvec| encrypt(&sk, &subvec))
        .collect::<Vec<_>>();

    (m, m_v, v, sk, x)
}

fn bench_whir_prover(c: &mut Criterion) {
    let mut rng = get_rng();

    let mut group = c.benchmark_group("whir_prover");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            // Number of variables is log2(D * height). We assume D*height is a power of two or
            // the workload is padded to the next power of two for the dummy MLE.
            let num_vars = (D * h).next_power_of_two().ilog2() as usize;
            let mle = MultilinearPolynomial::new(vec![F::from(1u64); 1 << num_vars], num_vars);

            let eval_point = (0..num_vars)
                .map(|i| F2::from(i as u64))
                .collect::<Vec<F2>>();

            let whir = Whir::<F2>::new(mle.num_variables(), &mut rng);

            let id = BenchmarkId::new("whir_prover", format!("T{}_H{}", t, h));
            group.bench_with_input(id, &(t, h), move |b, _inp| {
                // We measure just the prove call on the dummy mle
                b.iter(|| whir.prove(&mle, &eval_point));
            });
        }
    }
    group.finish();
}

fn bench_plaintext_matvec(c: &mut Criterion) {
    let mut group = c.benchmark_group("plaintext_matvec");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            let id = BenchmarkId::new("plaintext_matvec", format!("T{}_H{}", t, h));
            group.bench_with_input(id, &(t, h), |b, &(t, h)| {
                b.iter_batched(
                    || {
                        // Generate matrix data
                        let integer_width = D * t;
                        let integer_height = h;

                        let m_data = (1..=integer_height * integer_width)
                            .map(|x| x % 16)
                            .collect::<Vec<_>>();

                        // Generate vector data
                        let v = (1..=integer_width).map(|x| x % 16).collect::<Vec<_>>();
                        (m_data, v, integer_width)
                    },
                    |(m, v, w)| {
                        m.chunks(w)
                            .map(|row| {
                                row.iter()
                                    .zip(v.iter())
                                    .fold(0, |acc, (elem, rhs_elem)| acc + *elem * *rhs_elem)
                            })
                            .collect::<Vec<_>>()
                            .iter()
                            .map(|x| x % 16)
                            .collect::<Vec<_>>()
                    },
                    BatchSize::SmallInput,
                )
            });
        }
    }
    group.finish();
}

fn bench_vector_encryption(c: &mut Criterion) {
    let mut group = c.benchmark_group("vector_encryption");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            let id = BenchmarkId::new("vector_encryption", format!("T{}_H{}", t, h));
            group.bench_with_input(id, &(t, h), |b, &(t, h)| {
                b.iter_batched(
                    || {
                        let (_, _, v, sk, _) = setup_benchmark_data(t, h);
                        (v, sk)
                    },
                    |(v, sk)| {
                        v.par_chunks(D)
                            .map(|subvec| encrypt(&sk, &subvec))
                            .collect::<Vec<_>>()
                    },
                    BatchSize::SmallInput,
                )
            });
        }
    }
    group.finish();
}

fn bench_prover_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("prover_computation");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            let id = BenchmarkId::new("prover_computation", format!("T{}_H{}", t, h));

            // Preprocess once per parameter combination, then measure prove.
            let (m, _, _, _, x) = setup_benchmark_data(t, h);

            let (m_rq, m_polyring, m_mle, m_mle_over_base_f, z1_num_vars, mut transcript) =
                Prover::<D, F2>::preprocess(&m);

            group.bench_with_input(id, &(t, h), move |b, _| {
                b.iter(|| {
                    Prover::<D, F2>::prove(
                        &m_rq,
                        &m_polyring,
                        &m_mle,
                        &m_mle_over_base_f,
                        z1_num_vars,
                        &mut transcript,
                        &x,
                    )
                });
            });
        }
    }
    group.finish();
}

fn bench_result_decryption(c: &mut Criterion) {
    let mut group = c.benchmark_group("result_decryption");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            let id = BenchmarkId::new("result_decryption", format!("T{}_H{}", t, h));

            let (m, _, _v, sk, x) = setup_benchmark_data(t, h);

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

            group.bench_with_input(id, &(t, h), |b, _| {
                b.iter(|| {
                    proof
                        .y
                        .par_iter()
                        .map(|c| decrypt(&sk, c)[0])
                        .collect::<Vec<_>>()
                })
            });
        }
    }
    group.finish();
}

fn bench_verifier_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("verifier_computation");
    for &t in TS.iter() {
        for &h in HEIGHTS.iter() {
            let id = BenchmarkId::new("verifier_computation", format!("T{}_H{}", t, h));

            let (m, _, _v, _sk, x) = setup_benchmark_data(t, h);

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

            let field_elt_size = F2::field_size_in_bits() / 8;
            let sumcheck_1_size = proof.z1_sumcheck_proof.size_in_bytes();
            let sumcheck_2_size = proof.z3_sumcheck_proof.size_in_bytes();

            let proof_size = sumcheck_1_size
            + sumcheck_2_size
            + proof.r_mle_proof.proof.len()
            + field_elt_size // add one field elt for the claim
            + proof.m_mle_proof.proof.len()
            + field_elt_size; // add one field elt for the claim

            let (m_rq, z1_num_vars, transcript) = Verifier::<D, F2>::preprocess(&m);

            group.bench_with_input(id, &(t, h), |b, _| {
                b.iter(|| {
                    Verifier::<D, F2>::verify(
                        &m_rq,
                        z1_num_vars,
                        &mut transcript.clone(),
                        &x,
                        proof.clone(),
                    )
                })
            });
            println!("Proof size (in bytes): {}", proof_size);
        }
    }
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_plaintext_matvec,
        bench_vector_encryption,
        bench_prover_computation,
        bench_result_decryption,
        bench_verifier_computation,
        // bench_whir_prover,
);

criterion_main!(benches);
