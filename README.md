This repository provides an implementation of a verifiable plaintext matrix-private vector protocol.

__WARNING__ : This code was developped for academic purposes and is not ready for production use.

### Build
The build process mainly relies on Rust and has been tested on MacOS with ARM64 and Linux with AMD64.

#### Requirements

- [Rust](https://rust-lang.org/tools/install/)
- [Python](https://www.python.org/downloads/) (optional)

#### Instructions
To build the project, make sure Rust is properly installed and the repository is properly cloned.
Then, it suffices to run
```bash
cargo build --release
```
from the project's root directory.

### Tests and benchmarks
The project comes with unit and end-to-end tests as well as exhaustive benchmarks.
To execute the unit tests, run
```bash
cargo test --release
```

The end-to-end test can be executed with
```bash
cargo test --release -- --ignored
``` 

To reproduce the results presented in the paper, one may run the provided benchmarks by executing 
```bash
cargo bench --bench functionality_bench | tee <filename>
```
from the project root. This will save the benchmark results to the `<filename>` file.

The results can then be compiled and visualized using the provided Python script:
```bash
cd scripts
python3 gen_plots.py <filename>
```
This will generate `pgfplots` plots from the benchmarks.
The plots can be compiled by using your favourite LaTeX compiler to compile the [`main.tex`](scripts/main.tex) file.

### Code structure

- The [`arith`](src/arith) module defines the fields used in the protocol and provides arithmetic primitives related to (cyclotomic) polynomial rings and linear algebra.
- The [`rlwe`](src/rlwe) module provides basic functionalities related to Regev-style Ring-LWE encryption.
- The [`protocol`](src/protocol) module is divided into several submodules, as follows:
    - [`sumcheck`](src/protocol/sumcheck) provides an implementation of the sumcheck protocol for a product of linear sumchecks, inspired by the [Libra](https://doi.org/10.1007/978-3-030-26954-8_24) paper;
    - [`prover`](src/protocol/prover) and [`verifier`](src/protocol/verifier) respectively implement the SNARG protocol prover and verifier;
    - [`transcript`](src/protocol/transcript/) implements a Fiat-Shamir transcript using BLAKE3;
    - [`pcs`](src/protocol/pcs/) provides a wrapper around the [WHIR](https://github.com/WizardOfMenlo/whir) polynomial commitment scheme.
- The [`scripts`](scripts) folder provides Python utilities to interpret measures collected in benchmarks.

### Cryptographic parameters
Parameters related to the finite fields ($q$ extension degree) are defined in [`field.rs`](src/arith/field.rs).
Those specific to RLWE ciphertexts ($\sigma, \kappa$) are defined in the [`rlwe`](src/rlwe/mod.rs) module.
The rest ($N, m, n, p$) are defined directly in the [`functionality_bench.rs`](benches/functionality_bench.rs) file.

There is no need to change the parameters to reproduce the results presented in the paper. 
Executing the benchmark routine as described above will iterate through all relevant parameter sets.