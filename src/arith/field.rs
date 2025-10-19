use ark_ff::{Fp2, Fp2Config, Fp64, MontBackend, MontConfig, MontFp};

#[derive(MontConfig)]
#[modulus = "18446744069414584321"]
#[generator = "7"]
pub struct FConfig64;
pub type Field64 = Fp64<MontBackend<FConfig64, 1>>;

pub struct F2Config64;
impl Fp2Config for F2Config64 {
    type Fp = Field64;

    const NONRESIDUE: Self::Fp = MontFp!("7");

    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        MontFp!("1"),
        // Fq(7)**(((q^1) - 1) / 2)
        MontFp!("18446744069414584320"),
    ];
}

pub type Field64_2 = Fp2<F2Config64>;
