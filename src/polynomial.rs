use std::ops::Add;

const N: usize = 256;

pub trait PolynomialRing: Add + Sized {}

pub struct Polynomial {
    coeffs: Vec<i64>
}

impl From<Vec<i64>> for Polynomial{
    fn from(value: Vec<i64>) -> Self {
        Polynomial { coeffs: value }
    }
}

impl Polynomial{
    fn new
}

impl Add for &Polynomial {
    type Output = Polynomial;
    fn add(self, other: Self) -> Polynomial {
        let new_coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(a, b)| (a + b)) // à replacer par mod q
            .collect();
        Polynomial { coeffs: new_coeffs }
    }
}