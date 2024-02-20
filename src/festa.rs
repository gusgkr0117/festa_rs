//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use std::fmt;

use crate::{
    discrete_log::{bidlp, xgcd_big},
    ecFESTA::{evaluate_isogeny_chain, evaluate_isogeny_chain_for_basis, factored_kummer_isogeny},
    fields::{
        FpFESTA::Fp,
        FpFESTAExt::{
            BASIS_ORDER, D1, D1_FACTORED, D2, D2_FACTORED, DA, DA1, DA1_FACTORED, DA2,
            DA2_FACTORED, DA_FACTORED, L_POWER, M1, M2, THETA_ISOGENY_LENGTH, THETA_STRATEGY,
        },
    },
    pairing::{tate_pairing, weil_pairing},
    supersingular::{
        compute_canonical_kernel, entangled_torsion_basis, has_factored_order,
        isogeny_from_scalar_x_only, point_has_factored_order, random_isogeny_x_only, torsion_basis,
    },
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

pub struct TrapdoorInput {
    s1: BigUint,
    s2: BigUint,
    diagmatB: (BigUint, BigUint),
}

impl TrapdoorInput {
    pub fn new(s1: &BigUint, s2: &BigUint, diagmatB: &(BigUint, BigUint)) -> Self {
        TrapdoorInput {
            s1: s1.clone(),
            s2: s2.clone(),
            diagmatB: diagmatB.clone(),
        }
    }

    pub fn extract(&self) -> (BigUint, BigUint, (BigUint, BigUint)) {
        (self.s1.clone(), self.s2.clone(), self.diagmatB.clone())
    }
}

impl fmt::Display for TrapdoorInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({}, {}, ({}, {}))",
            self.s1, self.s2, self.diagmatB.0, self.diagmatB.1
        )
    }
}

pub struct FESTACrypto {
    E0: Curve,
    l_power: BigUint,
    d1: BigUint,
    d2: BigUint,
    d1d2: BigUint,
    clear_d1: BigUint,
    clear_d2: BigUint,
    clear_lb: BigUint,
    g2: BigUint,
    E0_lb_basis: (Point, Point),
    E0_d1_basis: (Point, Point),
}

impl FESTACrypto {
    pub fn new() -> Self {
        let E0 = Curve::new(&(Fq::TWO + Fq::FOUR));
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let basis_order = BigUint::from_slice(&BASIS_ORDER);
        let d1 = BigUint::from_slice(&D1);
        let d2 = BigUint::from_slice(&D2);
        let d1d2 = &d1 * &d2;
        let clear_d1 = xgcd_big(&d1, &l_power) * &d1;
        let clear_d2 = xgcd_big(&d2, &l_power) * &d2;
        let clear_lb = xgcd_big(&l_power, &d2) * &l_power;
        let E0_lb_basis = entangled_torsion_basis(&E0, &(basis_order / &l_power));
        let E0_d1_basis = torsion_basis(&E0, &D1_FACTORED, L_POWER as usize);
        let m1 = BigUint::from_slice(&M1);
        let m2 = BigUint::from_slice(&M2);
        let g1 = &m2 * BigUint::from_slice(&DA2) * &d2;
        let g1_inv = xgcd_big(&g1, &l_power);
        let g2 = &m1 * &d1 * g1_inv;

        FESTACrypto {
            E0,
            l_power,
            d1,
            d2,
            d1d2,
            clear_d1,
            clear_d2,
            clear_lb,
            g2,
            E0_lb_basis,
            E0_d1_basis,
        }
    }

    /// Generate a pair of public key and secret key of FESTA PKE
    pub fn keygen(&self) -> (PublicKey, SecretKey) {
        let mut rng = rand::thread_rng();
        // Setting the starting curve
        let E0 = self.E0;

        // Choose a random diagonal matrix A
        let (mut a11, mut a22) = (rng.gen_biguint(256), rng.gen_biguint(256));

        a11.set_bit(0, true);
        a22.set_bit(0, true);

        println!("[keygen] Pick random isogenies phiA1 and phiA2..");
        // Choose a random isogeny of degree dA1
        let (phiA1, EA_prime) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);

        // Choose a random isogeny of degree dA2
        let (phiA2, EA) = random_isogeny_x_only(&EA_prime, &DA2_FACTORED, 2);

        // Precompute some parameters
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let d1 = BigUint::from_slice(&D1);
        let d2 = BigUint::from_slice(&D2);
        let dA1 = BigUint::from_slice(&DA1);
        let dA2 = BigUint::from_slice(&DA2);
        let clear_d1 = xgcd_big(&d1, &l_power) * &d1;
        let clear_d2 = xgcd_big(&d2, &l_power) * &d2;
        let clear_lb = xgcd_big(&l_power, &d2) * &l_power;

        // Create torsion basis
        println!("[keygen] Create torsion basis of EA..");
        let (Pb, Qb) = self.E0_lb_basis;
        let (Pd1, Qd1) = self.E0_d1_basis;
        let (PA_prime_d2, QA_prime_d2) = torsion_basis(&EA_prime, &D2_FACTORED, L_POWER as usize);

        let (Pd1b, Qd1b) = (E0.add(&Pb, &Pd1), E0.add(&Qb, &Qd1));
        println!("[keygen] Evaluating the isogeny phiA1..");
        let (_, (imPd1b, imQd1b)) = evaluate_isogeny_chain_for_basis(
            &E0,
            &(Pd1b, Qd1b),
            &phiA1,
            &(&self.l_power * &self.d1),
            &dA1,
        );

        let (mut imPb, mut imQb) = (
            EA_prime.mul_big(&imPd1b, &clear_d1),
            EA_prime.mul_big(&imQd1b, &clear_d1),
        );

        let (Pbd2, Qbd2) = (
            EA_prime.add(&imPb, &PA_prime_d2),
            EA_prime.add(&imQb, &QA_prime_d2),
        );

        println!("[keygen] Evaluating the isogeny phiA2..");
        let (_, (imPbd2, imQbd2)) = evaluate_isogeny_chain_for_basis(
            &EA_prime,
            &(Pbd2, Qbd2),
            &phiA2,
            &(&self.l_power * &self.d2),
            &dA2,
        );

        (imPb, imQb) = (
            EA.mul_big(&imPbd2, &clear_d2),
            EA.mul_big(&imQbd2, &clear_d2),
        );

        let (imPAd2, imQAd2) = (
            EA.mul_big(&imPbd2, &clear_lb),
            EA.mul_big(&imQbd2, &clear_lb),
        );

        println!("[keygen] Generating the torsion basis of EA");
        let (PA_d2, QA_d2) = torsion_basis(&EA, &D2_FACTORED, L_POWER as usize);
        let ePQ = weil_pairing(&EA, &PA_d2, &QA_d2, &d2);
        println!("[keygen] Solving the BiDLP for the EA torsion basis");
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
        debug_assert!(
            point_has_factored_order(&EA, &R, &[(2, L_POWER)]),
            "R has wrong order"
        );
        debug_assert!(
            point_has_factored_order(&EA, &S, &[(2, L_POWER)]),
            "S has wrong order"
        );

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
    pub fn trapdoor_eval(&self, pk: &PublicKey, m: &TrapdoorInput) -> TrapdoorOutput {
        let (s1, s2, diagmatB) = m.extract();
        let (EA, R, S) = pk.extract();
        let (b11, b22) = diagmatB.clone();

        let (phi_1, E1) = isogeny_from_scalar_x_only(&self.E0, &D1_FACTORED, &s1, None);
        let (Pb, Qb) = self.E0_lb_basis;
        let (_, (imPb, imQb)) =
            evaluate_isogeny_chain_for_basis(&self.E0, &(Pb, Qb), &phi_1, &self.l_power, &self.d1);

        let (phi_2, E2) = isogeny_from_scalar_x_only(&EA, &D2_FACTORED, &s2, None);
        let (_, (imR, imS)) =
            evaluate_isogeny_chain_for_basis(&EA, &(R, S), &phi_2, &self.l_power, &self.d2);

        let (R1, S1) = (E1.mul_big(&imPb, &b11), E1.mul_big(&imQb, &b22));
        let (R2, S2) = (E2.mul_big(&imR, &b11), E2.mul_big(&imS, &b22));

        debug_assert!(
            point_has_factored_order(&E1, &R1, &[(2, L_POWER)]),
            "R1 has wrong order"
        );
        debug_assert!(
            point_has_factored_order(&E1, &S1, &[(2, L_POWER)]),
            "S1 has wrong order"
        );
        debug_assert!(
            point_has_factored_order(&E2, &R2, &[(2, L_POWER)]),
            "R2 has wrong order"
        );
        debug_assert!(
            point_has_factored_order(&E2, &S2, &[(2, L_POWER)]),
            "S2 has wrong order"
        );

        TrapdoorOutput::new(&E1, &R1, &S1, &E2, &R2, &S2)
    }

    /// Given the (2^b, 2^b)-isogeny between elliptic products and pairs
    /// of points L1_preimage, L2_preimage compute L1 and L2 by
    /// phi(L1_preimage) and Phi(L2_preimage)
    fn compute_Li_from_richelot_chain() -> (Point, Point) {
        todo!();
    }

    /// Wrapper function which recovers s1, s2 and the masking matrix B given
    /// the output of the (l,l)-chain during the inverse trapdoor function
    /// calculation
    fn recover_si_and_B(
        &self,
        E: &Curve,
        L1: &Point,
        L2: &Point,
        im_basis_bd1_E0: &(Point, Point),
        im_basis_d2_EA: &(Point, Point),
    ) -> TrapdoorInput {
        let imPd1d2 = E.mul_big(&L1, &self.l_power);
        let imQd1d2 = E.mul_big(&L2, &self.l_power);

        let imPb = E.mul_big(&L1, &self.d1d2);
        let (imPd1, imQd1) = (E.mul_big(&imPd1d2, &self.d2), E.mul_big(&imQd1d2, &self.d2));
        let (imPd2, imQd2) = (E.mul_big(&imPd1d2, &self.d1), E.mul_big(&imQd1d2, &self.d1));

        let (imPd1b, imQd1b) = im_basis_bd1_E0;

        let (phi_A1_Pb, phi_A1_Qb) = (
            E.mul_big(&imPd1b, &self.clear_d1),
            E.mul_big(&imQd1b, &self.clear_d1),
        );

        let (phi_A1_Pd1, phi_A1_Qd1) = (
            E.mul_big(&imPd1b, &self.l_power),
            E.mul_big(&imQd1b, &self.l_power),
        );

        let (phi_A2_dual_Pd2, phi_A2_dual_Qd2) = im_basis_d2_EA;

        let s1 = compute_canonical_kernel(E, &imPd1, &imQd1, &D1_FACTORED);
        let s2 = compute_canonical_kernel(E, &imPd2, &imQd2, &D2_FACTORED);
        println!("s1 : {}", s1);
        println!("s2 : {}", s2);

        todo!();
    }

    /// The FESTA inverted trapdoor function following algorithm 7 in the FESTA paper.
    /// Given a tuple c and the secret key, sk, this function recovers the integers
    /// s1, s2 and the masking matrix B
    pub fn trapdoor_inverse(&self, c: &TrapdoorOutput, sk: &SecretKey) -> TrapdoorInput {
        let (EA_prime, im_basis_bd1_E0, im_basis_d2_EA, A) = sk.extract();
        let (E1, R1, S1, E2, R2, S2) = c.extract();
        let (a11, a22) = A;
        let (R2_prime, S2_prime) = (
            E2.mul_big(&R2, &xgcd_big(&a11, &self.l_power)),
            E2.mul_big(&S2, &xgcd_big(&a22, &self.l_power)),
        );

        let (Pd1_1, Qd1_1) = torsion_basis(&E1, &D1_FACTORED, L_POWER as usize);
        let (Pd2_2, Qd2_2) = torsion_basis(&E2, &D2_FACTORED, L_POWER as usize);

        let (glue_P1, glue_Q1) = (R1, S1);
        let (glue_P2, glue_Q2) = (
            E2.mul_big(&R2_prime, &self.g2),
            E2.mul_big(&S2_prime, &self.g2),
        );

        let eP1Q1 = weil_pairing(
            &E1,
            &E1.double_iter(&glue_P1, 2),
            &E1.double_iter(&glue_Q1, 2),
            &BigUint::from(2u32).pow(THETA_ISOGENY_LENGTH),
        );
        let eP2Q2 = weil_pairing(
            &E2,
            &E2.double_iter(&glue_P2, 2),
            &E2.double_iter(&glue_Q2, 2),
            &BigUint::from(2u32).pow(THETA_ISOGENY_LENGTH),
        );
        debug_assert_eq!(
            eP1Q1 * eP2Q2,
            Fq::ONE,
            "The kerenl points are not maximally isotropic"
        );

        let (L1_1, L1_2) = (E1.add(&Pd1_1, &R1), Pd2_2);
        let (L2_1, L2_2) = (E1.add(&Qd1_1, &S1), Qd2_2);

        let P1P2 = CouplePoint::new(&glue_P1, &glue_P2);
        let Q1Q2 = CouplePoint::new(&glue_Q1, &glue_Q2);

        let L1 = CouplePoint::new(&L1_1, &L1_2);
        let L2 = CouplePoint::new(&L2_1, &L2_2);
        let image_points = [L1, L2];
        let E1E2 = EllipticProduct::new(&E1, &E2);

        println!("E1 : {}", E1.j_invariant());
        println!("E2 : {}", E2.j_invariant());

        let (E3E4, images) = product_isogeny(
            &E1E2,
            &P1P2,
            &Q1Q2,
            &image_points,
            THETA_ISOGENY_LENGTH as usize,
            &THETA_STRATEGY,
        );

        let (E3, E4) = E3E4.curves();
        println!("EA_prime : {}", EA_prime.j_invariant());
        println!("E3 : {}", E3.j_invariant());
        println!("E4 : {}", E4.j_invariant());
        let (P1, P2) = images[0].points();
        let (P3, P4) = images[1].points();

        let (L1, L2) = if E3.j_invariant() == EA_prime.j_invariant() {
            assert!(
                EA_prime.get_constant().equals(&E3.get_constant()) != 0,
                "{} + {} = {}",
                EA_prime.get_constant(),
                E3.get_constant(),
                EA_prime.get_constant() + E3.get_constant()
            );
            if (EA_prime.get_constant() + E3.get_constant()) == Fq::ZERO {
                ()
            }
            (P1, P3)
        } else {
            assert!(EA_prime.get_constant().equals(&E4.get_constant()) != 0);
            (P2, P4)
        };

        self.recover_si_and_B(&EA_prime, &L1, &L2, &im_basis_bd1_E0, &im_basis_d2_EA)
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use num_bigint::RandBigInt;
    use rand::thread_rng;

    #[test]
    fn run_festa_trapdoor() {
        let mut timelapse = Instant::now();
        let festa = FESTACrypto::new();
        println!(
            "festa initialization elapsed : {} ms",
            timelapse.elapsed().as_millis()
        );

        timelapse = Instant::now();
        let (pubkey, seckey) = festa.keygen();
        println!("keygen elapsed : {} ms", timelapse.elapsed().as_millis());

        // Pick random constants
        let d1 = BigUint::from_slice(&D1);
        let d2 = BigUint::from_slice(&D2);
        let mut rng = thread_rng();
        let s1 = rng.gen_biguint_range(&BigUint::from(0u32), &d1);
        let s2 = rng.gen_biguint_range(&BigUint::from(0u32), &d2);
        let (mut a11, mut a22) = (
            rng.gen_biguint(L_POWER as u64),
            rng.gen_biguint(L_POWER as u64),
        );
        a11.set_bit(0, true);
        a22.set_bit(0, true);
        let trapdoor_input = TrapdoorInput::new(&s1, &s2, &(a11, a22));

        println!("Input : {}", trapdoor_input);

        // Evaluating the FESTA trapdoor
        timelapse = Instant::now();
        let result = festa.trapdoor_eval(&pubkey, &trapdoor_input);
        println!(
            "trapdoor_eval elapsed : {} ms",
            timelapse.elapsed().as_millis()
        );

        // Inverting the FESTA trapdoor
        timelapse = Instant::now();
        let inverted = festa.trapdoor_inverse(&result, &seckey);
        println!(
            "trapdoor_inverse elapsed : {} ms",
            timelapse.elapsed().as_millis()
        );

        println!("Output : {}", inverted);
    }

    #[test]
    fn run_keygen() {
        let mut timelapse = Instant::now();
        let festa = FESTACrypto::new();
        println!(
            "festa initialization elapsed : {} ms",
            timelapse.elapsed().as_millis()
        );

        timelapse = Instant::now();
        let (pubkey, seckey) = festa.keygen();
        println!("keygen elapsed : {} ms", timelapse.elapsed().as_millis());
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
