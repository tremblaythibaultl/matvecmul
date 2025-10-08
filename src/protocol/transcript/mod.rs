use ark_crypto_primitives::sponge::{
    CryptographicSponge,
    poseidon::{PoseidonConfig, PoseidonSponge},
};
use ark_ff::{Field, PrimeField};

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
        self.sponge
            .absorb(&e.to_base_prime_field_elements().collect::<Vec<_>>());
    }

    pub fn squeeze(&mut self) -> F {
        let challenge: Vec<F::BasePrimeField> = self.sponge.squeeze_field_elements(1);

        self.sponge.absorb(&challenge);

        F::from_base_prime_field_elems(challenge).unwrap()
    }
}
