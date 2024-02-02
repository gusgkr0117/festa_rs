//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use crate::{
    fields::{
        FpFESTA::Fp,
        FpFESTAExt::{P0x, P0y, Q0x, Q0y, BASIS_ORDER, DA, DA1, DA2},
    },
    thetaFESTA::{
        product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, KummerLineIsogeny, Point,
    },
};

use num_bigint::{BigUint, RandBigInt};

/// Publickey structure of FESTA
pub struct PublicKey {
    EA: Curve,
    RA: Point,
    SA: Point,
}

/// Secretkey structure of FESTA
pub struct SecretKey {
    phiA: KummerLineIsogeny,
    diagmatA: (BigUint, BigUint),
}

/// Generate a pair of public key and secret key of FESTA PKE
pub fn keygen() {
    let mut rng = rand::thread_rng();
    // Setting the starting curve
    let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));

    // Choose a random diagonal matrix A
    let (a11, a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

    // Choose a random basis of order dA
    let basis_P0 = Point::new_xy(&Fq::new(&P0x, &Fp::ZERO), &Fq::new(&P0y, &Fp::ZERO));
    let basis_Q0 = Point::new_xy(&Fq::new(&Q0x, &Fp::ZERO), &Fq::new(&Q0y, &Fp::ZERO));
    
    let da = BigUint::from_slice(&DA);
    let basis_order = BigUint::from_slice(&BASIS_ORDER);
    let cofactor_da = basis_order/da;
    let basis_PA = E0.mul_big(&basis_P0, &cofactor_da);
    let mut basis_QA = E0.mul_big(&basis_Q0, &cofactor_da);

    let random_secret = rng.gen_biguint(256);
    basis_QA = E0.mul_big(&basis_QA, &random_secret);

    // Compute the kernel of order dA1
    let da1_kernel = E0.mul_big(&basis_QA, &BigUint::from_slice(&DA2));
}
