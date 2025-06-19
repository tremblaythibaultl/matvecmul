use ark_ff::PrimeField;

use crate::arith::cyclotomic_ring::CyclotomicRing;

/// Row-major matrix
pub struct Matrix<F: PrimeField> {
    data: Vec<F>,
    width: usize,
}

impl<F: PrimeField> Matrix<F> {
    /// Create a new row-major matrix from a vector of coefficients and a specified width.
    /// The length of the vector must be a multiple of the width.
    pub fn from_vec(data: Vec<F>, width: usize) -> Self {
        assert!(
            data.len() % width == 0,
            "Row-major representation length must be a multiple of width"
        );
        Matrix { data, width }
    }

    fn mul_vec(&self, rhs: &Vec<F>) -> Self {
        todo!()
    }

    /// Create a new matrix from a vector of coefficients and a specified width.
    /// The length of the vector must be a multiple of the width.
    pub fn process(&self) -> Self {
        let mut data_pr = Vec::<F>::with_capacity(self.data.len());

        for row in self.data.chunks(self.width) {
            data_pr.push(row[0]);
            for i in (1..self.width).rev() {
                data_pr.push(-row[i]);
            }
        }

        Matrix {
            data: data_pr,
            width: self.width,
        }
    }

    /// Convert the matrix to a vector of cyclotomic ring elements.
    /// Each row of the matrix is interpreted as a cyclotomic ring elements.
    /// The width of the matrix must match the degree of the cyclotomic polynomial.
    pub fn to_cyclotomic_relts<const D: usize>(&self) -> Vec<CyclotomicRing<D, F>> {
        let mut res = Vec::<CyclotomicRing<D, F>>::with_capacity(self.data.len() / self.width);

        for row in self.data.chunks(self.width) {
            res.push(CyclotomicRing::from_coeffs(&row));
        }

        res
    }
}

#[cfg(test)]
mod test {
    use crate::arith::{field::Field64, linalg::Matrix};

    type F = Field64;
    const HEIGHT: usize = 4;
    const WIDTH: usize = 4;

    fn pp(m: &Matrix<F>) {
        let mut iter = m.data.iter();
        for _ in 0..HEIGHT {
            for _ in 0..WIDTH {
                print!("{} ", iter.next().unwrap());
            }
            print! {"\n"}
        }
    }

    #[test]
    fn process() {
        // Use a simple, predictable matrix for testing
        let data = vec![
            F::from(1),
            F::from(2),
            F::from(3),
            F::from(4), // Row 1: [1, 2, 3, 4]
            F::from(5),
            F::from(6),
            F::from(7),
            F::from(8), // Row 2: [5, 6, 7, 8]
            F::from(9),
            F::from(10),
            F::from(11),
            F::from(12), // Row 3: [9, 10, 11, 12]
            F::from(13),
            F::from(14),
            F::from(15),
            F::from(16), // Row 4: [13, 14, 15, 16]
        ];

        let m = Matrix { data, width: WIDTH };

        let m_pr = m.process();

        // Expected result: each row should have first element unchanged,
        // then remaining elements in reverse order and negated
        let expected = vec![
            F::from(1),
            -F::from(4),
            -F::from(3),
            -F::from(2), // Row 1: [1, -4, -3, -2]
            F::from(5),
            -F::from(8),
            -F::from(7),
            -F::from(6), // Row 2: [5, -8, -7, -6]
            F::from(9),
            -F::from(12),
            -F::from(11),
            -F::from(10), // Row 3: [9, -12, -11, -10]
            F::from(13),
            -F::from(16),
            -F::from(15),
            -F::from(14), // Row 4: [13, -16, -15, -14]
        ];

        assert_eq!(m_pr.data, expected);
        assert_eq!(m_pr.width, WIDTH);
    }
}
