use std::marker::PhantomData;

use ark_crypto_primitives::sponge::{
    CryptographicSponge,
    poseidon::{PoseidonConfig, PoseidonSponge},
};
use ark_ff::{BigInteger, Field, PrimeField};
use ark_std::iterable::Iterable;
use sha3::{
    Shake256, Shake256Core,
    digest::{ExtendableOutput, Update, XofReader, core_api::CoreWrapper},
};

pub struct ShakeTranscript<F: Field> {
    hasher: CoreWrapper<Shake256Core>,
    _pd: PhantomData<F>,
}

impl<F: Field> ShakeTranscript<F> {
    pub fn new() -> Self {
        ShakeTranscript {
            hasher: Shake256::default(),
            _pd: PhantomData,
        }
    }

    pub fn absorb(&mut self, e: &F) {
        for elt in e.to_base_prime_field_elements() {
            let bytes = elt.into_bigint().to_bytes_le();
            self.hasher.update(&bytes)
        }
    }

    //TODO: need to squeeze_slice
    pub fn squeeze(&mut self, num: usize) -> Vec<F> {
        //TODO: Verify is this clone makes sense
        let mut reader = self.hasher.clone().finalize_xof();

        let ext_deg = F::extension_degree() as usize;

        // TODO: This assumes the field elements are 64 bits
        let mut output = vec![0u8; 8 * ext_deg * num];

        reader.read(&mut output);

        output
            .chunks(8 * ext_deg)
            .map(|bytes| {
                F::from_base_prime_field_elems(
                    bytes
                        .chunks(8)
                        .map(F::BasePrimeField::from_be_bytes_mod_order)
                        .collect::<Vec<F::BasePrimeField>>(),
                )
                .unwrap()
            })
            .collect::<Vec<F>>()
    }
}

pub struct PoseidonTranscript<F: Field> {
    sponge: PoseidonSponge<F::BasePrimeField>,
}

impl<F: Field> PoseidonTranscript<F>
where
    F: PrimeField + ark_crypto_primitives::sponge::Absorb,
{
    pub fn new(config: PoseidonConfig<F>) -> Self {
        Self {
            sponge: PoseidonSponge::new(&config),
        }
    }

    pub fn absorb(&mut self, e: &F) {
        self.sponge.absorb(&e);
    }

    pub fn squeeze(&mut self) -> F {
        let challenge = self.sponge.squeeze_field_elements(1);

        self.sponge.absorb(&challenge);

        challenge[0]
    }
}
