//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use std::time::Instant;

use crate::{
    discrete_log::{bidlp, xgcd_big},
    ecFESTA::{evaluate_isogeny_chain, factored_kummer_isogeny},
    fields::{
        FpFESTA::Fp,
        FpFESTAExt::{
            P0x, P0y, Q0x, Q0y, BASIS_ORDER, D1, D1_FACTORED, D2, D2_FACTORED, DA, DA1,
            DA1_FACTORED, DA2, DA2_FACTORED, DA_FACTORED, L_POWER,
        },
    },
    pairing::tate_pairing,
    supersingular::{isogeny_from_scalar_x_only, random_isogeny_x_only, torsion_basis},
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

    pub fn extract(&self) -> (Curve, Point, Point) {
        (self.EA, self.RA, self.SA)
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

    pub fn extract(&self) -> (Curve, (Point, Point), (Point, Point), (BigUint, BigUint)) {
        (
            self.EA_prime,
            self.im_basis_bd1_E0,
            self.im_basis_d2_EA,
            self.diagmatA.clone(),
        )
    }
}

pub struct TrapdoorOutput {
    E1: Curve,
    R1: Point,
    S1: Point,
    E2: Curve,
    R2: Point,
    S2: Point,
}

impl TrapdoorOutput {
    pub fn new(E1: &Curve, R1: &Point, S1: &Point, E2: &Curve, R2: &Point, S2: &Point) -> Self {
        TrapdoorOutput {
            E1: *E1,
            R1: *R1,
            S1: *S1,
            E2: *E2,
            R2: *R2,
            S2: *S2,
        }
    }

    pub fn extract(&self) -> (Curve, Point, Point, Curve, Point, Point) {
        (self.E1, self.R1, self.S1, self.E2, self.R2, self.S2)
    }
}

pub struct FESTACrypto {
    E0: Curve,
    l_power: BigUint,
    d1: BigUint,
    d2: BigUint,
    clear_d1: BigUint,
    clear_d2: BigUint,
    clear_lb: BigUint,
    E0_lb_basis: (Point, Point),
}

impl FESTACrypto {
    pub fn new() -> Self {
        let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let d1 = BigUint::from_slice(&D1);
        let d2 = BigUint::from_slice(&D2);
        let clear_d1 = xgcd_big(&d1, &l_power) * &d1;
        let clear_d2 = xgcd_big(&d2, &l_power) * &d2;
        let clear_lb = xgcd_big(&l_power, &d2) * &l_power;
        let E0_lb_basis = torsion_basis(&E0, &[(2u32, L_POWER)], 0);
        FESTACrypto {
            E0,
            l_power,
            d1,
            d2,
            clear_d1,
            clear_d2,
            clear_lb,
            E0_lb_basis,
        }
    }

    /// Generate a pair of public key and secret key of FESTA PKE
    pub fn keygen(&self) -> (PublicKey, SecretKey) {
        let mut rng = rand::thread_rng();
        // Setting the starting curve
        let E0 = self.E0;

        // Choose a random diagonal matrix A
        let (a11, a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

        // Choose a random isogeny of degree dA1
        let (phiA1, EA_prime) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);

        // Choose a random isogeny of degree dA2
        let (phiA2, EA) = random_isogeny_x_only(&EA_prime, &DA2_FACTORED, 2);

        // Precompute some parameters
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let d1 = BigUint::from_slice(&D1);
        let d2 = BigUint::from_slice(&D2);
        let clear_d1 = xgcd_big(&d1, &l_power) * &d1;
        let clear_d2 = xgcd_big(&d2, &l_power) * &d2;
        let clear_lb = xgcd_big(&l_power, &d2) * &l_power;

        // Create torsion basis
        let (Pb, Qb) = self.E0_lb_basis;
        let (Pd1, Qd1) = torsion_basis(&E0, &D1_FACTORED, L_POWER as usize);
        let (PA_prime_d2, QA_prime_d2) = torsion_basis(&EA_prime, &D2_FACTORED, L_POWER as usize);

        let (Pd1b, Qd1b) = (E0.add(&Pb, &Pd1), E0.add(&Qb, &Qd1));
        let ((_, imPd1b), (_, imQd1b)) = (
            evaluate_isogeny_chain(&E0, &Pd1b, &phiA1),
            evaluate_isogeny_chain(&E0, &Qd1b, &phiA1),
        );

        let (mut imPb, mut imQb) = (
            EA_prime.mul_big(&imPd1b, &clear_d1),
            EA_prime.mul_big(&imQd1b, &clear_d1),
        );

        let (Pbd2, Qbd2) = (
            EA_prime.add(&imPb, &PA_prime_d2),
            EA_prime.add(&imQb, &QA_prime_d2),
        );

        let ((_, imPbd2), (_, imQbd2)) = (
            evaluate_isogeny_chain(&EA_prime, &Pbd2, &phiA2),
            evaluate_isogeny_chain(&EA_prime, &Qbd2, &phiA2),
        );

        (imPb, imQb) = (
            EA.mul_big(&imPbd2, &clear_d2),
            EA.mul_big(&imQbd2, &clear_d2),
        );

        let (imPAd2, imQAd2) = (
            EA.mul_big(&imPbd2, &clear_lb),
            EA.mul_big(&imQbd2, &clear_lb),
        );

        let (PA_d2, QA_d2) = torsion_basis(&EA, &D2_FACTORED, L_POWER as usize);
        let ePQ = tate_pairing(&EA, &PA_d2, &QA_d2, &d2);
        let (a1, b1) = bidlp(&EA, &PA_d2, &imPAd2, &imQAd2, &D2_FACTORED, Some(ePQ));
        let (a2, b2) = bidlp(&EA, &QA_d2, &imPAd2, &imQAd2, &D2_FACTORED, Some(ePQ));

        let imPA_prime_d2 = EA_prime.add(
            &EA_prime.mul_big(&PA_prime_d2, &a1),
            &EA_prime.mul_big(&QA_prime_d2, &b1),
        );
        let imQA_prime_d2 = EA_prime.add(
            &EA_prime.mul_big(&PA_prime_d2, &a2),
            &EA_prime.mul_big(&QA_prime_d2, &b2),
        );

        let im_basis_bd1_E0 = (imPd1b, imQd1b);
        let im_basis_d2_EA = (imPA_prime_d2, imQA_prime_d2);

        // Mask torsion points using the matrix A
        let (R, S) = (EA.mul_big(&imPb, &a11), EA.mul_big(&imQb, &a22));

        // Set public key and secret key
        let public_key = PublicKey::new(&EA, &R, &S);
        let secret_key = SecretKey::new(&EA_prime, &im_basis_bd1_E0, &im_basis_d2_EA, &(a11, a22));

        (public_key, secret_key)
    }

    /// FESTA trapdoor function
    /// Input : the public key (EA, R, S)
    ///         two integers s1, s2
    ///         the masking matrix B
    /// Output : two elliptic curves E1 and E2 along with masked torsion bases
    ///         <R1, S1> = E1[l^b] and <R2, S2> = E2[l^b]
    pub fn trapdoor_eval(
        &self,
        pk: &PublicKey,
        s1: &BigUint,
        s2: &BigUint,
        diagmatB: &(BigUint, BigUint),
    ) -> TrapdoorOutput {
        let (EA, R, S) = pk.extract();
        let (b11, b22) = diagmatB.clone();

        let (phi_1, E1) = isogeny_from_scalar_x_only(&self.E0, &D1_FACTORED, s1, None);
        let (Pb, Qb) = self.E0_lb_basis;
        let ((_, imPb), (_, imQb)) = (
            evaluate_isogeny_chain(&self.E0, &Pb, &phi_1),
            evaluate_isogeny_chain(&self.E0, &Qb, &phi_1),
        );

        let (phi_2, E2) = isogeny_from_scalar_x_only(&EA, &D2_FACTORED, s2, None);
        let ((_, imR), (_, imS)) = (
            evaluate_isogeny_chain(&EA, &R, &phi_2),
            evaluate_isogeny_chain(&EA, &S, &phi_2),
        );

        let (R1, S1) = (E1.mul_big(&imPb, &b11), EA.mul_big(&imQb, &b22));
        let (R2, S2) = (E2.mul_big(&imR, &b11), EA.mul_big(&imS, &b22));

        TrapdoorOutput::new(&E1, &R1, &S1, &E2, &R2, &S2)
    }

    /// The FESTA inverted trapdoor function following algorithm 7 in the FESTA paper.
    /// Given a tuple c and the secret key, sk, this function recovers the integers
    /// s1, s2 and the masking matrix B
    pub fn trapdoor_inverse(
        &self,
        c: &TrapdoorOutput,
        sk: &SecretKey,
    ) -> (BigUint, BigUint, (BigUint, BigUint)) {
        let (EA_prime, im_basis_bd1_E0, im_basis_d2_EA, A) = sk.extract();
        let (E1, R1, S1, E2, R2, S2) = c.extract();
        let (a11, a22) = A;
        let (R2_prime, S2_prime) = (xgcd_big(&a11, &self.l_power), xgcd_big(&a22, &self.l_power));

        let (Pd1_1, Qd1_1) = torsion_basis(&E1, &D1_FACTORED, L_POWER as usize);
        let (Pd2_2, Qd2_2) = torsion_basis(&E2, &D2_FACTORED, L_POWER as usize);

        let (glue_P1, glue_Q1) = (R1, S1);
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_keygen() {
        let festa = FESTACrypto::new();
        festa.keygen();
    }

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
