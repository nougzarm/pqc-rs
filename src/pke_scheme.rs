use std::marker::PhantomData;

use crate::constants::PolyParams;
use crate::conversion::ByteEncode;
use crate::hash::{g, prf};
use crate::polynomial::{Polynomial, PolynomialNTT};

pub struct K_PKE<P: PolyParams> {
    k: usize,
    eta_1: usize,
    eta_2: usize,
    d_u: usize,
    d_v: usize,
    _marker: std::marker::PhantomData<P>,
}

impl<P: PolyParams> K_PKE<P> {
    pub fn new(k: usize, eta_1: usize, eta_2: usize, d_u: usize, d_v: usize) -> Self {
        K_PKE::<P> {
            k,
            eta_1,
            eta_2,
            d_u,
            d_v,
            _marker: PhantomData::<P>,
        }
    }

    /// Algorithm 13 : K-PKE.KeyGen(d)
    /// 
    /// Input : randomness d in B^32
    /// Output : (ek, dk) pair of encryption-decryption keys
    /// with : ek in B^(384*k + 32), and dk in B^(384*k)
    pub fn key_gen(&self, d: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
        let mut d_tmp = d.to_vec();
        d_tmp.extend_from_slice(&[self.k as u8]);
        let (rho, gamma) = g(&d_tmp);

        let mut n_var = 0usize;

        let a_ntt: Vec<Vec<PolynomialNTT<P>>> = vec![];
        for i in 0..self.k {
            let mut tmp_line = vec![];
            for j in 0..self.k {
                let mut input = [0u8; 34];
                input[0..32].copy_from_slice(&rho);
                input[32] = j as u8;
                input[33] = i as u8;
                tmp_line.push(PolynomialNTT::<P>::sample_ntt(&input));
            }
        };

        let mut s: Vec<Polynomial<P>> = vec![];
        for i in 0..self.k {
            s.push(Polynomial::<P>::sample_poly_cbd(&prf(self.eta_1, &gamma, &[n_var as u8]), self.eta_1));
            n_var += 1;
        };

        let mut e: Vec<Polynomial<P>> = vec![];
        for i in 0..self.k {
            e.push(Polynomial::<P>::sample_poly_cbd(&prf(self.eta_1, &gamma, &[n_var as u8]), self.eta_1));
            n_var += 1;
        };

        let s_ntt: Vec<PolynomialNTT<P>> = s.iter().map(|poly| poly.to_ntt()).collect();
        let e_ntt: Vec<PolynomialNTT<P>> = e.iter().map(|poly| poly.to_ntt()).collect();

        let mut t_ntt: Vec<PolynomialNTT<P>> = Vec::with_capacity(self.k);
        for i in 0..self.k {
            let mut pol_temp = PolynomialNTT::<P>::from(vec![0i64; P::N]);
            
            for j in 0..self.k {
                let product = &a_ntt[i][j] * &s_ntt[j];
                pol_temp = &pol_temp + &product; 
            }
            
            let t_i = &pol_temp + &e_ntt[i];
            t_ntt.push(t_i);
        }

        const CONST_D: usize = 12;

        let mut ek = Vec::new();
        for poly in &t_ntt {
            ek.extend(ByteEncode(&poly.coeffs, CONST_D));
        }
        ek.extend_from_slice(&rho);

        let mut dk = Vec::new();
        for poly in &s_ntt {
            dk.extend(ByteEncode(&poly.coeffs, CONST_D));
        }

        (ek, dk)
    }
}
