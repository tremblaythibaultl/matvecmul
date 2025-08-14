use criterion::{Criterion, criterion_group, criterion_main};
use matvec::arith::{
    cyclotomic_ring::CyclotomicRing, field::Field64, ntt::tfhe_based_ntt::TfheBasedNtt,
};

type Ring<const D: usize> = CyclotomicRing<D, Field64, TfheBasedNtt>;

fn bench_basic_mul<const D: usize>(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) {
    group.bench_function(format!("basic_mul {}", D), |b| {
        let x = Ring::<D>::get_random_bin();
        let y = Ring::<D>::get_random_bin();
        b.iter(|| {
            let _ = x.basic_mul(&y);
        });
    });
}

fn bench_tfhe_ntt_mul<const D: usize>(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) {
    group.bench_function(format!("tfhe ntt_mul {}", D), |b| {
        let x = Ring::<D>::get_random_bin();
        let y = Ring::<D>::get_random_bin();
        b.iter(|| {
            let _ = x.ntt_mul(&y);
        });
    });
}

macro_rules! for_sizes {
    ($group:expr, $f:ident, [$($d:expr),* $(,)?]) => {
        $( $f::<$d>($group); )*
    };
}

fn cyclotomic_ring_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("CyclotomicRing::mul");

    for_sizes!(&mut group, bench_basic_mul, [256, 512, 1024]);

    for_sizes!(&mut group, bench_tfhe_ntt_mul, [256, 512, 1024, 2048, 4096]);
}

criterion_group!(benches, cyclotomic_ring_mul);
criterion_main!(benches);
