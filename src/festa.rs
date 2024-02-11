//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use std::time::Instant;

use crate::{
    ecFESTA::{evaluate_isogeny_chain, factored_kummer_isogeny},
    fields::{
        FpFESTA::Fp,
        FpFESTAExt::{
            P0x, P0y, Q0x, Q0y, BASIS_ORDER, D1, D1_FACTORED, D2_FACTORED, DA, DA1, DA1_FACTORED,
            DA2, DA2_FACTORED, DA_FACTORED, L_POWER,
        },
    },
    supersingular::{random_isogeny_x_only, torsion_basis},
    thetaFESTA::{
        product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, KummerLineIsogeny, Point,
    },
};

use num_bigint::{BigUint, RandBigInt};
use rand::random;

/// Publickey structure of FESTA
pub struct PublicKey {
    EA: Curve,
    RA: Point,
    SA: Point,
}

impl PublicKey {
    pub fn new(EA: &Curve, RA: &Point, SA: &Point) -> Self {
        PublicKey {
            EA: *EA,
            RA: *RA,
            SA: *SA,
        }
    }
}

/// Secretkey structure of FESTA
pub struct SecretKey {
    EA_prime: Curve,
    im_basis_bd1_E0: (Point, Point),
    im_basis_d2_EA: (Point, Point),
    diagmatA: (BigUint, BigUint),
}

impl SecretKey {
    pub fn new(
        EA_prime: &Curve,
        im_basis_bd1_E0: &(Point, Point),
        im_basis_d2_EA: &(Point, Point),
        diagmatA: &(BigUint, BigUint),
    ) -> Self {
        SecretKey {
            EA_prime: *EA_prime,
            im_basis_bd1_E0: *im_basis_bd1_E0,
            im_basis_d2_EA: *im_basis_d2_EA,
            diagmatA: diagmatA.clone(),
        }
    }
}

/// Generate a pair of public key and secret key of FESTA PKE
pub fn keygen() {
    let mut rng = rand::thread_rng();
    // Setting the starting curve
    let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));

    // Choose a random diagonal matrix A
    let (a11, a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

    // Choose a random isogeny of degree dA1
    let isogeny_start_time = Instant::now();
    let (phiA1, EA_prime) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);
    println!(
        "random_isogeny_x_only elapsed : {} ms",
        isogeny_start_time.elapsed().as_millis()
    );
    // Choose a random isogeny of degree dA2
    let (phiA2, EA) = random_isogeny_x_only(&EA_prime, &DA2_FACTORED, 2);

    // Create a 2^b-torsion basis
    let l_power = BigUint::from(2u32).pow(L_POWER);
    let d1 = BigUint::from_slice(&D1);
    let (Pb, Qb) = torsion_basis(&E0, &[(2u32, L_POWER)], 0);
    let (Pd1, Qd1) = torsion_basis(&E0, &D1_FACTORED, L_POWER as usize);
    let (PA_prime_d2, QA_prime_d2) = torsion_basis(&EA_prime, &D2_FACTORED, L_POWER as usize);

    let (Pd1b, Qd1b) = (E0.add(&Pb, &Pd1), E0.add(&Qb, &Qd1));
    let ((_, imPd1b), (_, imQd1b)) = (
        evaluate_isogeny_chain(&E0, &Pd1b, &phiA1),
        evaluate_isogeny_chain(&E0, &Qd1b, &phiA1),
    );

    let (imPb, imQb) = (
        EA_prime.mul_big(&imPd1b, &d1),
        EA_prime.mul_big(&imQd1b, &d1),
    );

    let (Pbd2, Qbd2) = (
        EA_prime.add(&imPb, &PA_prime_d2),
        EA_prime.add(&imQb, &QA_prime_d2),
    );

    let ((_, imPbd2), (_, imQbd2)) = (
        evaluate_isogeny_chain(&EA_prime, &Pbd2, &phiA2),
        evaluate_isogeny_chain(&EA_prime, &Qbd2, &phiA2),
    );

    let (imPAd2, imQAd2) = (
        EA_prime.mul_big(&imPbd2, &l_power),
        EA_prime.mul_big(&imQbd2, &l_power),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keygen_torsion_basis() {
        let mut rng = rand::thread_rng();
        // Setting the starting curve
        let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));

        // Choose a random diagonal matrix A
        let (a11, a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

        // Choose a random isogeny of degree dA1
        let mut timewatch = Instant::now();
        let (phiA1, EA_prime) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);
        println!(
            "generate random da1 elapsed : {} ms",
            timewatch.elapsed().as_millis()
        );

        timewatch = Instant::now();
        // Choose a random isogeny of degree dA2
        let (phiA2, EA) = random_isogeny_x_only(&EA_prime, &DA2_FACTORED, 2);
        println!(
            "generate random da2 elapsed : {} ms",
            timewatch.elapsed().as_millis()
        );

        // Create a 2^b-torsion basis
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let d1 = BigUint::from_slice(&D1);
        timewatch = Instant::now();
        let (Pb, Qb) = torsion_basis(&E0, &[(2u32, L_POWER)], 0);
        println!(
            "generate l^e torsion basis of E0 elapsed : {} ms",
            timewatch.elapsed().as_millis()
        );

        timewatch = Instant::now();
        let (Pd1, Qd1) = torsion_basis(&E0, &D1_FACTORED, L_POWER as usize);
        println!(
            "generate d1 torsion basis of E0 elapsed : {} ms",
            timewatch.elapsed().as_millis()
        );

        timewatch = Instant::now();
        let (PA_prime_d2, QA_prime_d2) = torsion_basis(&EA_prime, &D2_FACTORED, L_POWER as usize);
        println!(
            "generate d2 torsion basis of EA elapsed : {} ms",
            timewatch.elapsed().as_millis()
        );
    }
}
