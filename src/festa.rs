//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use crate::{
    ecFESTA::factored_kummer_isogeny, fields::{
        FpFESTA::Fp,
        FpFESTAExt::{P0x, P0y, Q0x, Q0y, BASIS_ORDER, D1, DA, DA1, DA1_FACTORED, DA2, DA2_FACTORED, DA_FACTORED, L_POWER},
    }, supersingular::{random_isogeny_x_only, torsion_basis}, thetaFESTA::{
        product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, KummerLineIsogeny, Point,
    }
};

use num_bigint::{BigUint, RandBigInt};
use rand::random;

/// Publickey structure of FESTA
pub struct PublicKey {
    EA: Curve,
    RA: Point,
    SA: Point,
}

/// Secretkey structure of FESTA
pub struct SecretKey {
    phiA: Vec<KummerLineIsogeny>,
    diagmatA: (BigUint, BigUint),
}



/// Generate a pair of public key and secret key of FESTA PKE
pub fn keygen() {
    let mut rng = rand::thread_rng();
    // Setting the starting curve
    let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));

    // Choose a random diagonal matrix A
    let (a11, a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

    // Choose a random isogeny of degree dA1
    let (phi1, EA1) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);

    // Choose a random isogeny of degree dA2
    let (phi2, EA2) = random_isogeny_x_only(&EA1, &DA2_FACTORED, 2);

    // Create a 2^b-torsion basis
    let l_power = BigUint::from(2u32).pow(L_POWER);
    let d1 = BigUint::from_slice(&D1);
    let (Pb, Qb) = torsion_basis(&E0, &l_power, 0);
    let (Pd1, Qd1) = torsion_basis(&E0, &d1, L_POWER as usize);

    let (Pd1b, Qd1b) = (E0.add(&Pb, &Pd1), E0.add(&Qb, &Qd1));
    
}

#[cfg(test)]
mod tests{
    use super::keygen;

    #[test]
    fn run_keygen() {
        keygen();
    }
}