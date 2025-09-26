use ark_ff::{Field, PrimeField};

use crate::{
    arith::{cyclotomic_ring::CyclotomicRing, polynomial_ring::PolynomialRing, ring::Ring},
    protocol::sumcheck::multilinear::MultilinearPolynomial,
    rlwe::RLWE,
};

/// Row-major matrix over a ring R
#[derive(Debug)]
pub struct Matrix<R: Ring> {
    data: Vec<R>,
    width: usize,
}

impl<R: Ring> Matrix<R> {
    /// Create a new row-major matrix from a vector of ring elements and a specified width.
    /// The length of the vector must be a multiple of the width.
    pub fn from_vec(data: Vec<R>, width: usize) -> Self {
        assert!(
            data.len() % width == 0,
            "Row-major representation length must be a multiple of width"
        );
        Matrix { data, width }
    }

    fn mul_vec(&self, _rhs: &Vec<R>) -> Self {
        todo!()
    }

    /// Rearrange the columns of the matrix such that when rows are interpreted as elements of the cyclotomic ring of degree `d`,
    /// the constant coefficient of the inner product of ring elements yields the same result as the inner product of the original elements.
    /// The length of the vector must be a multiple of the width.
    pub fn process(&self, d: usize) -> Self {
        let mut data_pr = Vec::<R>::with_capacity(self.data.len());

        for row in self.data.chunks(d) {
            data_pr.push(row[0].clone());
            for i in (1..d).rev() {
                data_pr.push(row[i].neg());
            }
        }

        Matrix {
            data: data_pr,
            width: self.width,
        }
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.data.len() / self.width
    }
}

impl<F: PrimeField> Matrix<F> {
    /// Lift the matrix over `F` a prime field to a matrix over `Rq` a cyclotomic ring of degree `D`.
    /// The columns are rearranged to satisfy the relevant inner product property.
    pub fn lift_to_rq<const D: usize>(&self) -> Matrix<CyclotomicRing<D, F>> {
        assert!(
            self.data.len() % self.width == 0,
            "Matrix length must be a multiple of width"
        );

        assert!(
            self.width % D == 0,
            "Matrix width must be a multiple of the degree of the cyclotomic ring"
        );

        // Rearrange the columns of the matrix
        let m_pr = self.process(D);

        let data = m_pr
            .data
            .chunks(D)
            .map(|row| CyclotomicRing::<D, F>::from_coeffs(row))
            .collect();

        Matrix {
            data,
            width: self.width / D,
        }
    }

    /// Multiplies a matrix over a prime field by a vector over the same field.
    pub fn mat_vec_mul(&self, rhs: &Vec<F>) -> Vec<F> {
        assert!(
            self.width == rhs.len(),
            "Matrix width must match vector length"
        );
        let result = self
            .data
            .chunks(self.width)
            .map(|row| {
                row.iter()
                    .zip(rhs.iter())
                    .fold(F::ZERO, |acc, (elem, rhs_elem)| acc + *elem * *rhs_elem)
            })
            .collect::<Vec<_>>();
        result
    }
}

impl<const D: usize, F: PrimeField> Matrix<CyclotomicRing<D, F>> {
    pub fn lift_to_polynomial_ring(&self) -> Matrix<PolynomialRing<D, F>> {
        let data = self
            .data
            .iter()
            .map(|elem| PolynomialRing::<D, F>::from_cyclotomic(elem))
            .collect();

        Matrix {
            data,
            width: self.width,
        }
    }

    pub fn to_mle_evals(&self) -> (Vec<F>, usize) {
        let evals = self
            .data
            .iter()
            .map(|elem| elem.coeffs.clone())
            .flatten()
            .collect::<Vec<_>>();

        let num_variables = (self.data.len() * D).next_power_of_two().ilog2() as usize;

        (evals, num_variables)
    }
}

impl<R: Ring> Matrix<R> {
    /// Multiplies a matrix over a ring by a vector over of RLWE ciphertexts over the same ring.
    pub fn mat_rlwe_vec_mul(&self, rhs: &Vec<RLWE<R>>) -> Vec<RLWE<R>> {
        assert!(
            self.width == rhs.len(),
            "Matrix width must match vector length"
        );

        let result = self
            .data
            .chunks(self.width)
            .map(|row| {
                row.iter()
                    .zip(rhs.iter())
                    .fold(RLWE::<R>::zero(), |mut acc, (elem, rhs_elem)| {
                        acc.add_assign(&rhs_elem.mul_constant(&elem));
                        acc
                    })
            })
            .collect::<Vec<_>>();

        result
    }
}

// TODO: replace this with properly traited implementation.
pub fn gaussian_elimination<F: Field>(matrix: &mut [Vec<F>]) -> Vec<F> {
    let size = matrix.len();
    assert_eq!(size, matrix[0].len() - 1);

    for i in 0..size - 1 {
        for j in i..size - 1 {
            echelon(matrix, i, j);
        }
    }

    for i in (1..size).rev() {
        eliminate(matrix, i);
    }

    for i in 0..size {
        if matrix[i][i] == F::zero() {
            println!("Infinitely many solutions");
        }
    }

    let mut result: Vec<F> = vec![F::zero(); size];
    for i in 0..size {
        result[i] = matrix[i][size] / matrix[i][i];
    }
    result
}

fn echelon<F: Field>(matrix: &mut [Vec<F>], i: usize, j: usize) {
    let size = matrix.len();
    if matrix[i][i] == F::zero() {
    } else {
        let factor = matrix[j + 1][i] / matrix[i][i];
        (i..size + 1).for_each(|k| {
            let tmp = matrix[i][k];
            matrix[j + 1][k] -= factor * tmp;
        });
    }
}

fn eliminate<F: Field>(matrix: &mut [Vec<F>], i: usize) {
    let size = matrix.len();
    if matrix[i][i] == F::zero() {
    } else {
        for j in (1..i + 1).rev() {
            let factor = matrix[j - 1][i] / matrix[i][i];
            for k in (0..size + 1).rev() {
                let tmp = matrix[i][k];
                matrix[j - 1][k] -= factor * tmp;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use ark_ff::AdditiveGroup;

    use super::Matrix;
    use crate::arith::{cyclotomic_ring::CyclotomicRing, field::Field64};

    type F = Field64;
    type TestRing = CyclotomicRing<4, F>;
    const HEIGHT: usize = 4;
    const WIDTH: usize = 4;

    fn pp(m: &Matrix<TestRing>) {
        let mut iter = m.data.iter();
        for _ in 0..HEIGHT {
            for _ in 0..WIDTH {
                print!("{:?} ", iter.next().unwrap().coeffs);
            }
            print! {"\n"}
        }
    }

    #[test]
    fn process() {
        // Create test ring elements
        let data = vec![
            TestRing::from_coeffs(&[F::from(1), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(2), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(3), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(4), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(5), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(6), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(7), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(8), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(9), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(10), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(11), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(12), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(13), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(14), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(15), F::ZERO, F::ZERO, F::ZERO]),
            TestRing::from_coeffs(&[F::from(16), F::ZERO, F::ZERO, F::ZERO]),
        ];

        let m = Matrix { data, width: WIDTH };
        let m_pr = m.process(4);

        // Verify that the first element of each row is unchanged
        // and the remaining elements are negated
        assert_eq!(m_pr.data[0].coeffs[0], F::from(1));
        assert_eq!(m_pr.data[1].coeffs[0], -F::from(4));
        assert_eq!(m_pr.data[4].coeffs[0], F::from(5));
        assert_eq!(m_pr.data[5].coeffs[0], -F::from(8));

        assert_eq!(m_pr.width, WIDTH);
    }
}
