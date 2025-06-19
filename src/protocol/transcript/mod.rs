use std::marker::PhantomData;

use ark_ff::{BigInteger, Field, PrimeField};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

#[derive(Clone)]
pub struct Blake3Transcript<F: Field> {
    hasher: blake3::Hasher,
    _pd: PhantomData<F>,
}

impl<F: Field> Blake3Transcript<F> {
    pub fn new() -> Self {
        Blake3Transcript {
            hasher: blake3::Hasher::new(),
            _pd: PhantomData,
        }
    }

    pub fn absorb(&mut self, e: &F) {
        for elt in e.to_base_prime_field_elements() {
            let bytes = elt.into_bigint().to_bytes_le();
            self.hasher.update(&bytes);
        }
    }

    pub fn absorb_bytes(&mut self, bytes: &[u8]) {
        self.hasher.update(bytes);
    }

    pub fn absorb_bytes_par(&mut self, bytes: &[u8]) {
        self.hasher.update_rayon(bytes);
    }

    // not sure if this makes sense - find a way to hash in parallel for better verifier times
    pub fn absorb_rayon(&mut self, slice: &[F]) {
        let bytes: Vec<u8> = slice
            .par_iter()
            .flat_map(|elt| {
                elt.to_base_prime_field_elements()
                    .map(|base_prime_elt| base_prime_elt.into_bigint().to_bytes_le())
                    .collect::<Vec<Vec<u8>>>()
            })
            .flatten()
            .collect();

        // only makes sense if bytes is > 128 KiB
        self.hasher.update_rayon(&bytes);
    }

    pub fn squeeze(&mut self, num: usize) -> Vec<F> {
        let mut output_reader = self.hasher.finalize_xof();

        let ext_deg = F::extension_degree() as usize;

        // TODO: This assumes the field elements are 64 bits
        let mut output = vec![0u8; 8 * ext_deg * num];

        output_reader.fill(&mut output);
        self.absorb_bytes(&output);
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
