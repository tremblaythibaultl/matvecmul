use std::{
    collections::{HashMap, hash_map::Entry},
    sync::{Arc, OnceLock, RwLock},
};

use ark_ff::PrimeField;
use tfhe_ntt::prime64::Plan;

use crate::arith::ntt::ntt::Ntt;

// Inspired by https://github.com/zama-ai/tfhe-rs/blob/5fa8cc8563cfd87fba4ad1a0421df648c7cfeb82/tfhe/src/core_crypto/commons/math/ntt/ntt64.rs
// Key is (D, modulus).
type PlanMap = RwLock<HashMap<(usize, u64), Arc<OnceLock<Arc<Plan>>>>>;
static PLANS: OnceLock<PlanMap> = OnceLock::new();
fn plans() -> &'static PlanMap {
    PLANS.get_or_init(|| RwLock::new(HashMap::new()))
}

#[derive(Clone, Debug, Default)]
pub struct TfheBasedNtt;

impl TfheBasedNtt {
    // Gets a static plan for the 64 bit NTT. If the polynomial size or the prime modulus are not suitable, it panics.
    fn ntt_plan_64<F: PrimeField>(polynomial_size: usize) -> Arc<Plan> {
        assert_eq!(F::MODULUS.as_ref().len(), 1, "Modulus must fit in 64 bits.");
        let modulus = F::MODULUS.as_ref()[0];
        let global_plans = plans();
        let get_plan = || {
            let plans = global_plans.read().unwrap();
            let plan = plans.get(&(polynomial_size, modulus)).cloned();
            drop(plans);
            plan.map(|p| {
                p.get_or_init(|| {
                    Arc::new(
                        Plan::try_new(polynomial_size, modulus)
                            .expect("Failed to create tfhe-ntt plan"),
                    )
                })
                .clone()
            })
        };

        if let Some(plan) = get_plan() {
            return plan;
        }

        // Could not find a plan of the given size, we lock the map again and try to insert it.
        let mut plans = global_plans.write().unwrap();
        if let Entry::Vacant(v) = plans.entry((polynomial_size, modulus)) {
            v.insert(Arc::new(OnceLock::new()));
        }

        drop(plans);
        get_plan().unwrap()
    }
}

impl Ntt for TfheBasedNtt {
    fn mul<F: PrimeField>(lhs: &[F], rhs: &[F]) -> Vec<F> {
        assert_eq!(
            lhs.len(),
            rhs.len(),
            "NTT multiplication requires equal sized inputs."
        );
        let polynomial_size = lhs.len();
        let plan = Self::ntt_plan_64::<F>(polynomial_size);
        // From https://github.com/zama-ai/tfhe-rs/blob/main/tfhe-ntt/examples/mul_poly_prime.rs
        // Converting to u64 should be safe as we assume that the modulus fits in 64 bits and that's checked at NTT plan creation.
        let mut lhs_u64 = lhs
            .iter()
            .map(|c| c.into_bigint().as_ref()[0])
            .collect::<Vec<_>>();
        let mut rhs_u64 = rhs
            .iter()
            .map(|c| c.into_bigint().as_ref()[0])
            .collect::<Vec<_>>();
        plan.fwd(&mut lhs_u64);
        plan.fwd(&mut rhs_u64);
        plan.mul_assign_normalize(&mut lhs_u64, &rhs_u64);
        plan.inv(&mut lhs_u64);
        lhs_u64.iter().map(|&x| F::from(x)).collect::<Vec<_>>()
    }
}
