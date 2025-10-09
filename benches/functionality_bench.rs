use criterion::{Criterion, criterion_group, criterion_main};
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

pub const D: usize = 1024;
pub const INTEGER_WIDTH: usize = 2 * D;
pub const INTEGER_HEIGHT: usize = 2 * D;
pub type F = Field64;
pub type F2 = Field64_2;

fn setup_benchmark_data() -> (
    Matrix<F>,
    Vec<F>,
    Vec<u64>,
    Vec<CyclotomicRing<D, F>>,
    Vec<RLWE<CyclotomicRing<D, F>>>,
) {
    // Generate matrix data
    let data = (1..=INTEGER_HEIGHT * INTEGER_WIDTH)
        .map(|x| F::from(x as u64))
        .collect::<Vec<_>>();
    let m = Matrix::from_vec(data, INTEGER_WIDTH);

    // Generate vector data
    let v = (1..=INTEGER_WIDTH)
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
    // Create a dummy MLE of size equal to the size of `r_mle` in the protocol
    let num_vars = (D * INTEGER_HEIGHT).ilog2() as usize;
    let mle = MultilinearPolynomial::new(vec![F::from(1u64); 1 << num_vars], num_vars);

    let eval_point = (0..num_vars)
        .map(|i| F2::from(i as u64))
        .collect::<Vec<F2>>();

    let mut rng = get_rng();
    let whir = Whir::<F2>::new(mle.num_variables(), &mut rng);

    c.bench_function("whir_prover", |b| b.iter(|| whir.prove(&mle, &eval_point)));
}

fn bench_matrix_creation(c: &mut Criterion) {
    c.bench_function("matrix_creation", |b| {
        b.iter(|| {
            let data = (1..=INTEGER_HEIGHT * INTEGER_WIDTH)
                .map(|x| F::from(x as u64))
                .collect::<Vec<_>>();
            Matrix::from_vec(data, INTEGER_WIDTH)
        });
    });
}

fn bench_plaintext_matvec(c: &mut Criterion) {
    let (m, _, v, _, _) = setup_benchmark_data();
    let v_field = v.iter().map(|&x| F::from(x)).collect::<Vec<F>>();

    c.bench_function("plaintext_matvec", |b| {
        b.iter(|| m.mat_vec_mul(&v_field));
    });
}

fn bench_vector_encryption(c: &mut Criterion) {
    let (_, _, v, sk, _) = setup_benchmark_data();

    c.bench_function("vector_encryption", |b| {
        b.iter(|| {
            v.chunks(D)
                .map(|subvec| encrypt(&sk, &subvec))
                .collect::<Vec<_>>()
        });
    });
}

fn bench_prover_computation(c: &mut Criterion) {
    let (m, _, _, _, x) = setup_benchmark_data();

    let (m_rq, m_polyring, m_mle, z1_num_vars, mut transcript) = Prover::<D, F2>::preprocess(&m);

    c.bench_function("prover_computation", |b| {
        b.iter(|| {
            Prover::<D, F2>::prove(&m_rq, &m_polyring, &m_mle, z1_num_vars, &mut transcript, &x)
        });
    });
}

fn bench_result_decryption(c: &mut Criterion) {
    let (m, _, _v, sk, x) = setup_benchmark_data();

    let (m_rq, m_polyring, m_mle, z1_num_vars, mut transcript) = Prover::<D, F2>::preprocess(&m);

    let proof =
        Prover::<D, F2>::prove(&m_rq, &m_polyring, &m_mle, z1_num_vars, &mut transcript, &x);

    // only consider the mask for now
    c.bench_function("result_decryption", |b| {
        b.iter(|| {
            proof
                .y
                .iter()
                .map(|c| decrypt(&sk, c)[0])
                .collect::<Vec<_>>()
        });
    });
}

fn bench_verifier_computation(c: &mut Criterion) {
    let (m, _, _v, _sk, x) = setup_benchmark_data();

    let (m_rq, m_polyring, m_mle, z1_num_vars, mut transcript) = Prover::<D, F2>::preprocess(&m);

    let proof =
        Prover::<D, F2>::prove(&m_rq, &m_polyring, &m_mle, z1_num_vars, &mut transcript, &x);

    c.bench_function("verifier_computation", |b| {
        b.iter(|| Verifier::<D, F2>::verify(&m, &x, proof.clone()));
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_matrix_creation,
        bench_plaintext_matvec,
        bench_vector_encryption,
        bench_prover_computation,
        bench_result_decryption,
        bench_verifier_computation,
        bench_whir_prover
);

criterion_main!(benches);
