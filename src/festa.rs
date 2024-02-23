//! keygen, encrypt and decrypt of FESTA Public key encryption scheme
//! It's an implementation of FESTA-128
use std::fmt;

use crate::{
    discrete_log::{bidlp, xgcd_big},
    ecFESTA::{
        evaluate_isogeny_chain, evaluate_isogeny_chain_for_basis, factored_kummer_isogeny,
        CurveIsomorphism,
    },
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

    pub fn encode(&self) -> Vec<u8> {
        let EA = self.EA_prime.get_constant();
        let ((x1, y1), (x2, y2)) = (
            self.im_basis_bd1_E0.0.to_xy(),
            self.im_basis_bd1_E0.1.to_xy(),
        );
        let ((x3, y3), (x4, y4)) = (self.im_basis_d2_EA.0.to_xy(), self.im_basis_d2_EA.1.to_xy());
        let mut result: Vec<u8> = Vec::new();
        result.extend_from_slice(&EA.encode());
        result.extend_from_slice(&x1.encode());
        result.extend_from_slice(&y1.encode());
        result.extend_from_slice(&x2.encode());
        result.extend_from_slice(&y2.encode());
        result.extend_from_slice(&x3.encode());
        result.extend_from_slice(&y3.encode());
        result.extend_from_slice(&x4.encode());
        result.extend_from_slice(&y4.encode());
        let a11 = self.diagmatA.0.to_bytes_le();
        result.extend(&a11);
        result.extend(&vec![0; 82 - a11.len()]);
        let a22 = self.diagmatA.1.to_bytes_le();
        result.extend(&a22);
        result.extend(&vec![0; 82 - a22.len()]);
        result
    }

    pub fn decode(bytearray: &Vec<u8>) -> Self {
        let total_len = Fq::ENCODED_LENGTH * 9 + 82 * 2;
        let mut offset: usize = 0;
        let mut var = [Fq::ZERO; 9];
        debug_assert_eq!(bytearray.len(), total_len);

        for i in 0..9 {
            let r: u32;
            (var[i], r) = Fq::decode(&bytearray[offset..offset + Fq::ENCODED_LENGTH]);
            debug_assert_eq!(r, u32::MAX);
            offset += Fq::ENCODED_LENGTH;
        }
        let a11 = BigUint::from_bytes_le(&bytearray[offset..offset + 82]);
        offset += 82;
        let a22 = BigUint::from_bytes_le(&bytearray[offset..offset + 82]);
        Self {
            EA_prime: Curve::new(&var[0]),
            im_basis_bd1_E0: (
                Point::new_xy(&var[1], &var[2]),
                Point::new_xy(&var[3], &var[4]),
            ),
            im_basis_d2_EA: (
                Point::new_xy(&var[5], &var[6]),
                Point::new_xy(&var[7], &var[8]),
            ),
            diagmatA: (a11, a22),
        }
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

    pub fn encode(&self) -> Vec<u8> {
        let A1 = self.E1.get_constant();
        let (R1x, R1y) = self.R1.to_xy();
        let (S1x, S1y) = self.S1.to_xy();
        let A2 = self.E2.get_constant();
        let (R2x, R2y) = self.R2.to_xy();
        let (S2x, S2y) = self.S2.to_xy();
        let mut result = Vec::new();
        result.extend_from_slice(&A1.encode());
        result.extend_from_slice(&R1x.encode());
        result.extend_from_slice(&R1y.encode());
        result.extend_from_slice(&S1x.encode());
        result.extend_from_slice(&S1y.encode());
        result.extend_from_slice(&A2.encode());
        result.extend_from_slice(&R2x.encode());
        result.extend_from_slice(&R2y.encode());
        result.extend_from_slice(&S2x.encode());
        result.extend_from_slice(&S2y.encode());
        result
    }

    pub fn decode(binary: &Vec<u8>) -> Self {
        let mut offset = 0;
        let mut var = [Fq::ZERO; 10];
        debug_assert_eq!(binary.len(), Fq::ENCODED_LENGTH * 10);

        for i in 0..10 {
            let r;
            (var[i], r) = Fq::decode(&binary[offset..offset + Fq::ENCODED_LENGTH]);
            debug_assert_eq!(r, u32::MAX);
            offset += Fq::ENCODED_LENGTH;
        }

        Self {
            E1: Curve::new(&var[0]),
            R1: Point::new_xy(&var[1], &var[2]),
            S1: Point::new_xy(&var[3], &var[4]),
            E2: Curve::new(&var[5]),
            R2: Point::new_xy(&var[6], &var[7]),
            S2: Point::new_xy(&var[8], &var[9]),
        }
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
    inv_scalar: BigUint,
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
        let inv_scalar = xgcd_big(&(&d1 * &d1d2 * &m1), &l_power);
        let g1 = &m2 * BigUint::from_slice(&DA2) * &d2;
        let g1_inv = xgcd_big(&g1, &l_power);
        let g2 = &m1 * &d1 * g1_inv;

        FESTACrypto {
            E0,
            l_power,
            d1,
            d2,
            d1d2,
            inv_scalar,
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
        let (mut a11, mut a22) = (
            rng.gen_biguint(L_POWER as u64),
            rng.gen_biguint(L_POWER as u64),
        );

        a11.set_bit(0, true);
        a22.set_bit(0, true);

        println!("[keygen] Pick random isogenies phiA1 and phiA2..");
        // Choose a random isogeny of degree dA1
        let (phiA1, EA_prime) = random_isogeny_x_only(&E0, &DA1_FACTORED, 2);

        // Choose a random isogeny of degree dA2
        let (phiA2, EA) = random_isogeny_x_only(&EA_prime, &DA2_FACTORED, 2);

        // Precompute some parameters
        let dA1 = BigUint::from_slice(&DA1);
        let dA2 = BigUint::from_slice(&DA2);

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
            EA_prime.mul_big(&imPd1b, &self.clear_d1),
            EA_prime.mul_big(&imQd1b, &self.clear_d1),
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
            EA.mul_big(&imPbd2, &self.clear_d2),
            EA.mul_big(&imQbd2, &self.clear_d2),
        );

        let (imPAd2, imQAd2) = (
            EA.mul_big(&imPbd2, &self.clear_lb),
            EA.mul_big(&imQbd2, &self.clear_lb),
        );

        println!("[keygen] Generating the torsion basis of EA");
        let (PA_d2, QA_d2) = torsion_basis(&EA, &D2_FACTORED, L_POWER as usize);
        let ePQ = weil_pairing(&EA, &PA_d2, &QA_d2, &self.d2);
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

    /// Helper function to recover the masking matrix B from the
    /// output of the (l,l)-chain during the inverse trapdoor function
    /// calculation
    pub fn recover_diagB(
        &self,
        E: &Curve,
        phi_1_dual_R_scaled: &Point,
        phi_A1_Pb: &Point,
        phi_A1_Qb: &Point,
    ) -> (BigUint, BigUint) {
        let (debug_value, mut alpha) = bidlp(
            E,
            phi_1_dual_R_scaled,
            phi_A1_Pb,
            phi_A1_Qb,
            &[(2, L_POWER)],
            None,
        );

        assert!(
            phi_1_dual_R_scaled.equals(
                &(E.add(
                    &E.mul_big(phi_A1_Pb, &debug_value),
                    &E.mul_big(phi_A1_Qb, &alpha)
                ))
            ) != 0
        );
        println!("R = {}P + {}Q", debug_value, alpha);

        alpha = &self.inv_scalar * alpha;
        let beta = xgcd_big(&alpha, &self.l_power);

        (alpha, beta)
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

        let s1 = compute_canonical_kernel(
            E,
            &imPd1,
            &imQd1,
            &D1_FACTORED,
            Some((phi_A1_Pd1, phi_A1_Qd1)),
        );
        let s2 = compute_canonical_kernel(
            E,
            &imPd2,
            &imQd2,
            &D2_FACTORED,
            Some((*phi_A2_dual_Pd2, *phi_A2_dual_Qd2)),
        );
        println!("s1 : {}", s1);
        println!("s2 : {}", s2);

        let diagmatB = self.recover_diagB(E, &imPb, &phi_A1_Pb, &phi_A1_Qb);

        TrapdoorInput::new(&s1, &s2, &diagmatB)
    }

    /// The FESTA inverted trapdoor function following algorithm 7 in the FESTA paper.
    /// Given a tuple c and the secret key, sk, this function recovers the integers
    /// s1, s2 and the masking matrix B
    pub fn trapdoor_inverse(&self, c: &TrapdoorOutput, sk: &SecretKey) -> TrapdoorInput {
        let (EA_prime, im_basis_bd1_E0, im_basis_d2_EA, A) = sk.extract();
        let (E1, R1, S1, E2, R2, S2) = c.extract();
        let (a11, a22) = A;

        println!("E1 : {}", E1.j_invariant());
        println!("E2 : {}", E2.j_invariant());

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
            if E3.get_constant() == EA_prime.get_constant() {
                (P1, P3)
            } else {
                let isom = CurveIsomorphism::new(&E3, &EA_prime);
                (isom.eval(&P1), isom.eval(&P3))
            }
        } else {
            debug_assert!(EA_prime.j_invariant().equals(&E4.j_invariant()) != 0);
            if E4.get_constant() == EA_prime.get_constant() {
                (P2, P4)
            } else {
                let isom = CurveIsomorphism::new(&E4, &EA_prime);
                (isom.eval(&P2), isom.eval(&P4))
            }
        };

        self.recover_si_and_B(&EA_prime, &L1, &L2, &im_basis_bd1_E0, &im_basis_d2_EA)
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use crate::ecFESTA::CurveIsomorphism;

    use super::*;
    use num_bigint::RandBigInt;
    use rand::thread_rng;

    #[test]
    fn isomorphism_test() {
        let A1_str = "64b8e83c0af74fc401fd47b4136314a752e23263c7b9a2028427c384b3e2799fb87837d7e28d8eb2d226eb496afe34cd5dc1c043b6b1837a115f9d4e0745430288dd0331157ae4f94d7e9028f0db0e121b3353254cd75a61cd4ef2061a798ab86eea37e4958c58b8b7021eac47aaafe5f534bbdf0ac19b36adfa20e68ceec09a4819bfebd9ec38a308b64ef9eafced4515710467858fcdc17e71dc8f167ba416171300a41683195bada9e84927eda6cea790d01f40c556355b8234f1f7533ba3013d268f6aac60348ea7ca2c1b3378e8a76b7b8f71a908e92c3d2e71346b18b28196e5bed344df21ee47ea52ce937d51773c260d8ebd7ae18984ff451e762be800fc1a030df10d7b8f0edec670ed768af2d016b5e1be6ce08c3d9f44ac8faefd0f81d86206bd3d57214be1030b33be8bd8d138ce897af5b776dc7c74e3714b444bd26d08";
        let A2_str = "c9a46f8732873fa1b4c86326ce7ede73e256d4e03f098ca1a0d6023d6c86cb0f40a8f8cd655129996aa48854acba1d81f3a8e3bf4271c6728088dbdd138a3d9513979ad2a0eb5139a68378409f835289ebc57de2ee76d0287700269eef26b916f35e3510259c705ca55f788aab445ae29adeeee1a17f168db09bd38ce7ebc2bf5d65b9f6e4e2cf1ca16e7dfa2cec751e74ae362adbd025329bae47a487899f16780b087976a075c323980253e68a64def6ad1735efced056ff6f6c45ef9f56e087c7ee59827b3addcf8332d3dc17fd91bc08ea0f70289b74d307c9385b4f6fe4fbbdff0ad412d01d65a436455d9127dae7fe72688e75a7074d11a45a06ea8e6e1fdae22792b2ba087794c4177b4b5a7d7c9ff15c239d867659b68b558a1cdddf6d191b61c58096e3a108122c2c26f56687fc495f875e4f791df4b4bf659edecc82305604";
        let (A1, _) = Fq::decode(&hex::decode(A1_str).unwrap());
        let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());
        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);

        let order = BigUint::from(2u32).pow(L_POWER);
        let cofactor = BigUint::from_slice(&BASIS_ORDER) / &order;
        let (P, Q) = entangled_torsion_basis(&E1, &cofactor);

        let isom = CurveIsomorphism::new(&E1, &E2);
        let (imP, imQ) = (isom.eval(&P), isom.eval(&Q));

        assert!(E2.on_curve(&imP) && E2.on_curve(&imQ));
        assert!(E2.mul_big(&imP, &order).isinfinity() != 0);
        assert!(E2.mul_big(&imQ, &order).isinfinity() != 0);

        assert!(point_has_factored_order(&E2, &imP, &[(2u32, L_POWER)]));
        assert!(point_has_factored_order(&E2, &imQ, &[(2u32, L_POWER)]));
    }

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
        let (mut b11, mut b22) = (
            rng.gen_biguint(L_POWER as u64),
            rng.gen_biguint(L_POWER as u64),
        );
        b11.set_bit(0, true);
        b22.set_bit(0, true);
        let trapdoor_input = TrapdoorInput::new(&s1, &s2, &(b11, b22));

        println!("Input : {}", trapdoor_input);

        // Evaluating the FESTA trapdoor
        timelapse = Instant::now();
        let result = festa.trapdoor_eval(&pubkey, &trapdoor_input);
        println!(
            "trapdoor_eval elapsed : {} ms",
            timelapse.elapsed().as_millis()
        );

        println!(
            "trapdoor output : {}",
            result
                .encode()
                .iter()
                .map(|b| format!("{:02X}", b))
                .collect::<Vec<String>>()
                .join("")
        );
        println!(
            "secret key : {}",
            seckey
                .encode()
                .iter()
                .map(|b| format!("{:02X}", b))
                .collect::<Vec<String>>()
                .join("")
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

    fn hex_to_bytes(s: &str) -> Option<Vec<u8>> {
        if s.len() % 2 == 0 {
            (0..s.len())
                .step_by(2)
                .map(|i| {
                    s.get(i..i + 2)
                        .and_then(|sub| u8::from_str_radix(sub, 16).ok())
                })
                .collect()
        } else {
            None
        }
    }

    #[test]
    fn test_festa_trapdoor_inverse() {
        let secret_key_binary = "375041DAE5FCBF7B3C5FCB9724049159D1761B679B8A378C54165FF928AA26E66F223C35E1FFECBB6088CDC412DF46BE734F32660811C3D55823F91CADD47352EA1AC7F3C8E9614C12871D69276004EF5ED1E5D7A513860BE64DD775B567B7610D765B44A8397BC61D0B20CE71FD1D89643BDC39F0B05D417CB458EE2D4C75EC989ACA7518B4D866C93B67BD36CD96A4D2850111A26EE957F9AF39502F767EEBF90D919A7379F3883EEA1F5794E3BE6438C635F35238DBC44A7FD97516D253040821E659002FD43F0A851039E983223B7DA0DEFDA2568633D47C4D8234F8C98E46E6341E0C8CA7C81EE72B82AE49144BD99AD4B5BCB0F3A27CC5A67DBB78CC92F2B563F29DA9AA69FF55BF2FAD166C31C23BF7D14B85BC7AF3EB44BC1681FEF37C5FA2650BDC9EC3F169388D0BF839480502E5DBF2CCACD6E278DC947F7A3471055E35146B1567CC524BE5C79C72A164D49ABEBCB6F63E72B0F8B1E4C118F6C951437C44E9D625A2CF87D004BB678CC79760CEA35625B7CA2296FAAF702B6CB4C0986C78E1A48198542506623FBD18F1C7DCC5B4A1A09DC5F99BA56478A88DC6E1E60F116B34BD10560A5571DE0647A91DECD6BB9F8861AF0CD44069D5BC21AF2D21C6908DB678C75C5D527D1BFC38F4779FEE996F83EFFE988BE154E40C0F762EC040A58B12E00352A0875FEFB8D34EEF004641A0551FB6745D3CBA0DE59885CDC9389CBB14BF35E814D957F8B2A1D628F8CBE65AEF8B64E1407E46CD5A35C415B223B909021E241FA0B48F65CDACA523162071276D5FA7991E9F08594A891FDEC6B9EF16A8A82BA8BAC0C7C9E2FCFCB6D6FD981A08F77C092217628532B8143017DDEC4E3074546FF1FE0EA698AAFB4A57EC1B0E4DAEF790E3690E86CCF0C4DEDDD403CE9FE605E0091D76E19EDDA962E8625F41D5608B97E1504BD0140BC38FB390103DF6D049A213E6CF9C215BD79B2AA3AC6A1A4446B1A42900EE76C3D08856F317C2116C520FAFB42005AA990A1038C1AA0011320E86A648E7D463A14582D7736550D90599DEDC5A4529FB48E0D1DD31CAD9C9EE1B842FDC12F24B16AD4FA054C31C4878207FF512117AC62B4F09C074D142C31BA39D2E9BF99F954469B5C5A50611D4159C1002F1B34FB93555F1C41906757D661779053BE6ABC5E6DE653AE49A8683E26BD0CAE488C37E0ACF245FC8CA20F8FDF6FA519073D55BB67F64DA61A995AF459D1CF1EA601B83991517585CE39403B6DCD685227E1B7B121705B803D510EFA134844AE300644E061B8A6C9594F2CCB39DDC6CD8BC699BC7DDB11A9302DD46D0B2E7A7CF60EF9728957675AF0949AB86D0B6249F8FDF893A3CAA4C9149F9B6ADC8155E110D260637F9AFB62EC6A1EDC30DE6025090E190D21C84759B8B95A5F77E086E424AD256136FB3AE9D42AFBC56C1C121750638AC87A49ABB32E47101BE2DB92FDF95FC88599560BB4BCEBB02ED0A529223CC5028E7FF2D1D7542F5943A5DAA54A9054E5E9C6133DACD906D6F8B547CBA55909D1AE0246FE2E6D1549D0029A68677BA29D8DC99A7D35EA6586D7A021747786969FB6727CA8C5D239051FFF357CC635BDC036749C5136316FFD01ADE76620BFF38010351DDCD37D26C5D45917E8A3BA9A72E72EFD3138843589327552687C4FD07135531BB0701DFA5B1D181E1E644CF4B9D73C9818E79F94BA84EF742DBF79CB2BD97E7941888C999233339DD354F55C76A71F48EBED085E9A58C2A05AA983FF211E2B35D159357BB3D93924DE34403BADFE78E8973E0E05C93C9B0A30D9F79E42018E79413C144D47A1DC297CF22581DDE9C126FCFE88DCA6CDF134F97CC0EFC45368F5D25284D47FE4C7BADA5DFC085D8D9DC3A28702EC9D8709DCD1140F4FF985A304290694CED2C248E06E37A3C6740ECC67A2A3D3BAB8DFD28C40C1D07B0E3C639A3EE0180EB1FECA106EA8B544D08951326B9C494043F085E631617F9BA8E69A0A2E651A6AE31D06B17FDAF783DC88363E89F3AA63B1FF006659DCCBADF005EF19CFF9C8F5E19D4E57834635026376BC3AB05104878AE0454743729A292326A57A94B08393ADA4F2C89FEFECDA78809E4E606B4831F11D300C5B01AEAD9F4CA9E92827C37EAEB312AB2BA94E23D632A8CDD59C5A4DACC15752B9DD9E69BBE2910B321B9201FC39CEE02BD1DC351BB89B2D25C66AA1144FBD7663B6C5D20D7879F8B67F23F78D8216A0EAF7D41C633E67B3FAAEB75F84A47D78267D94E0EB1A6AC67574AC68F9B4D3F820278C1D83F29AEC0B83B0EB2D837E22C7326A9071304E78CA489A0AB8A8C4ABA941AA099CDFD513D872D60FE42B5DE4001CB1E6288BEFAFDA0207438058076B18F809C1449D8AADD83B1F3C16F9B7B775245DC6EF43B61034E54BDD8D5EDA2F92FA8AFE0DDC2DC8B98991719A9BC0F92CE6FAEA3B431E323F95E4830567E4CF326E978FA92F31680E8899F486031C002E4B8FEE5C4868EAD49B6FC312F4C145CE4796DF65FE5A8A79FF1D264140676C82FD48D6F6BFD8A8877449E3AA19897144C4B3E44A4850FE5FA8A9ECA7CD3EC33C28AB7DC45DC3D693BCE76627C53AC17B74A507C0735CE690CBC4964CBAC41C6EF74469A1494A8BE1362471240DF0682993C29ADC2B406BEA28501339AAFDC7253304F46FDC8DEAC7CB17293ED6663D5F66EFA8E7C1A2309D8EC0F97B5D29392315167DE18D104CA1582FF4B3B1D8459236C8166BC4819CB062AE1A70A1013045548386287338EA99624DC0B71FAB52CB640107A821E5AF0A60B82718466F8171E7A526CCD7DEEA5AA7EA74E8456070F593E1687E39EC48E910F0E9D7C818ECC62557D9301CB2B8C937251B8F2DB0B2C717B5FC4155F35DB8BF7C68CD2E0727723D5B4C9EE69BA1B365F14DD2F9EAAE1BDBBDC23E3E5896F5278267B687580640282FAD9E786696AC1A43A87AB21B7E9047F8EE4AB31484F1F0808A97FDAE49E0B9EA9ACBEEACF1D65A6E1C0FC49BED4EFA558FC4E755B538FE305C138CF28058AD1B7C5663E70C0FCBBD4B7B472CEB3017B370A6F3BD5898856EC112233C83A050EC36897880068C8057946397CE1E7E338A2B095EFECC9B14FC3C9A92136B9E63C61ECF14070D77A8BE93273E81072E0ECA0C793B58B6EC0D1FF65B88A4032AE729583E9D4E83929EF86F7FB2EF2243E29ED62498541501F828DD77C688881509B115BFF787C18B8AFD778A6491D053A850AB26ACED88BE7962F7F4E5E721AC0DE0D4AD343E5EFDD1EF98061A5182EC9AD6EEFD368F315CBAD4080DA1DDC34F0EF38C90BA71A9DE05AA1B99A7F82F5C949A023029389F0E0668B2D88B450FE3930FDEBADF8003C1737115C5318F1208820103FDD8BB679F5CC597611762CFFDF82D5C9F648440937FC8F5FFAC9B9F2A26BEDAAEEE1F594D4C621637D1DD92F98E5A76B0BBABB5977DA5EE8574D17C17C6B12B8296EEAFF888275F3678EB83C3ED217CE5E1F82B4879F06F101BDAD0EA46586CA7F9360E48707D38EBA50A938FF6B2BA4023195711139F24E94F1C627810437FA68EE068269B6785D09AAD6C913735A0016A5A95BB60FA7D65EF8F0D8D23C8E98872F782C45991DEFCF0B8441E5033C52354CE6E559D6BBC4E093AB316BC7A44E96E371BDD6458CE82D4556E0A5871696E2FC1312AB70D3E3D63D095A859B4CE9A78E73A9B93FA3C76E03F7B484769677F28E273ACA821105BC0E969FA7993708E862C196DBD08A62D7614C42C8EC0A82C4ECB620BF30C04FFCC4FAEEDD8D7940AE31A791D42A84A6A0C913F74CBAAC7A80E169333BD055BADEB6F703A2C5CEE341027D02CB74DFFA2A0314E2AB13FCF4221FD124F36E19DC2703ED5162D025676971E9332FB96A73DF698C251574D9D7E86B698510FF47344003FD8456EC312181B97458951156EDB142CEE42ADA1D7A16B76A2A100D8B6A03DFFDF3A0B4D014FD05AB25CF2213B68A8B1D042DD602E2B6FD989F54E353B97BE436741EC025C9D8C4C1AB5662BDCE405DF2B996691543E3EBB10E5D8B0E64DBA2342E48E6877783101B1C40976345FAD9FD327E39EE61A311718FEEB9579FA509A4FAD9F91BE87F5660D1FB5DEC230A1A4E2DBA51984A122697CCA4FF16CFB7497B659916355343871ACA11955825A57AEFE2DF277BBD717F1BDFF2930389B04FD8FD7CCF7401F8C0676095657D179FE48F17D78A4F3045B9026F428C3C69FAACFC813CE4BB986DD9BAE744CE00000001287E71851B866FD540C932F91DF70784371D4C277DDBB68FE83E35D15B2B1FF6CFE0D0E2AE2243BF70A54C5724509A767E81401E8AB12C74E7724F032042EFE6EC0663C69D05FB91824FA2BAFA79010000";
        let trapdoor_output_binary = "F646038ADF39416D7F1161BA6A568734EB07BD3E008CFDB35D45A793B69EE76E69837CCBFC4E93A5D30228B273D17C7045B992F822C7909C233ABA4BCB0A1049FF21077F87CD3DF11BC024C7C1DE4730AEBBBD62674693C3949A344BB7452CB509F0C6CCE29707745BDDC0037F893D9659462BCC0504B03C1882A144C7223ABADB567D5E31C574046CAB0A75FA8723D1B77B151E9CFF2C8969DFBBE142872B2B1C10CE7E27A8B80A836FDF4C5664E1985DBAF885D5501FF71147FDA3B90814B343F87B2D21AB52AE7BDDEB5109A323EC50996B5B2A8A3106B9AC68C3F56ACFDE7D98D44155B039F9D0D5BEFA88E8B1C421F8E971F53AB4573A4324F054A7171AB89CDC8E2A029C7C0979D229557CD67A15AAF78F5F13B71F37D5D20E18D9F1A4A770C2F16A11C8374FDFF965C68F8D76331FD3238BEDE75C805BFDDDDB58B002C03BFB1447938B6175DAB295BC8E1207DDE3DCE563C147EF46FB4AB9217ADF58B90B4E7EE19820169B64CE8CF92616861AA79DA5350DB18988F52D0A83D9BD7CD43667A4CBA1D66A86C24E96EB8EEE73405F4CC4F98728F6F5A909FDB7576A1F3FE5E312AEFC05EA1C1D8F2A9E52D4B139AB40914D890B8D15561BAF4C50946F53763C2F81FB3EFDC7247BEAF7111DB069BD5F7471FA22DE156FD167374AB8F7867D3E40790C7C4CF7C0688DEC5A1FB3162947986828E69382844261D8717324B8034A51CA9D20D083E516AD75E05314BBD3649A6D3344266EB517095AB61F5743C3D124F08A766D45B0F1D1410E9AEC2DED80B1C069CF3F679F6D053BF7A38DF8D46217F00D03802EEAFA68839EB569D5A893846D6078C19B73B5D826A83C34CC96FD5C6DB18B68AEAE75230952DA6273C8236758EB02C204DF7025BBADFE9C7AD278E3D38769133C0AAC47B9F4A9F51D319359BEBDF045BC3B6665E630C4DA4406326FDD92F873B2A71F959134F52BC0ACA80B346AE840E856A389D11C43549C4B6C9B337C974AC6CA7FE3C4DCE400F7FF4B16C6394F86D3F2B3554904F269A1B8C0E226ACE467987EAB9BAA48160D27EAAA9A840CAB17CC75823140397749E06537565F1B72904C3468DC79C2A1D3AC1A01CF5C235959BF4F0FF0CBB5721875FA39FD2A762701D4097444C8C36908215DE16A75286A415A4A945F0C5A3F511CF3101331CBD077F5341DFB6F3812DC0DA453C40614F89E466550536037372517CD71C2F8E1F75382BF9D030975713F7DA966B731BCC31F0C464C071829FF942A498C30B776504AF0DE418FB768F132CA9C19C103BA9AD460307D26CB80622FB3CBAC7EB09F42C574582185B04ED5B8075EE69DB6AFBA71A941095E0F90AA9A22AD331D33F70B86058B9915EA6F9F303BC4418BD6F1264B54D4AA68369151C1B9C12996D819F04D993672BA94AABB1D84AF6614BD6B5EC5A63DA94468106FB6F21D73A474302EFC8E8D1852B1E2E90B31B081F1CBCE71A95AAB400F0243946C8F5605A40D08362EA95C07D58E128B3B94D07AAE99610469B9B652D2E5B427E55E4D8DF71A0ADE15D992A666AE2F413DB64E8AC29975D023BA475BB62CE5A6ED8624E80E2E32BF4629C9614961135634454492ECA51977C616327A030ED5C941E0DD7AC898EA6088DD80D4CF85B90FDF0A6078D1FDD592D778BB49D78B4FC8F22042AC8CFC3A08665220F8E37974AA67C6CE46A120A8BB8EB90494A68E1B549AA89693F5B00CA07CE7568A0773D1FE9D0FF205896B9B1A5CEED8A0F1D6DBF545572B0A89EFB1AD71C9DEA390686E7360483361E364D502B54AD06E6CECA04C131AB151483AB2D561D9E4794511569C16B7B1721D13DB05B4D73878F8A8F8D862C786682FCE46FB87C9BFC2F506F4A3CA0EC69E9E2AA01A34FD75AD40863E3E7F0CB8AD7FE93F589FD22DD2FD16F22BA13F536316A69875DB8FC0474F92A11D8241D36839C95263BA9FE201EE3FB4AD60484003FD3ABCE016C2CB0AE753066B21DC8345AB0001631DB36752B5E38337B71FC8D8C562E5C666C402BA35D0A88B1F6D402C9794258D3AA308C8D76960C4BA4F006FC3A195601C75BBEBEEE9BF2DF3AE312D951ED0B21444940291EEA5CDFCA27AE6FF2849CF865E6EEB328BCA66BA7EDD6D33257057C5C128D95F3BCEFF825F1A0C64223C6BCBCC001A5527866A017369BD0BB60CC4C239B5AA380422F0776CAFCC989C81926384499F57B6154852956ED7A67FABB8944BB04EFDA39A33C009247494554F055EDD2C226267A5BC5A5A14F973883E4B44429617FCCB140AC8E791BB0D1A0CBDD2CCFE25695CDC9AAE63C7C9AC917124EAA173B69EBBD3D9FD8FC05CCDC20A217C96F0139B75A2B0A3DFEC87D7FE1D8EEE653881EFA21D98F13AD48B1F94A6F7A8079AC50382D3FF7512C4F00158DAC00A64F424D9FEF64511144FE1F3FCB9763C9A5DA60A2382A06ACA9D0CA2EE935EA28FE089F45330149B10B1E93B11582BBD0B22891121B12389B579C2318B991225497D3C9A70B9D339254C3C5F231666EE9BD51E99220ED3BA2DE6E0276AADEF29796BF0B81263E5D9D130743651C60208D164353353727CB452F2B09209080FFE0091681920DEC804CFFCD2B903C7D8F12DC65885ADF3AC0ECA418B85A52B6629440BF8179767880B67F348B52E3F9870756DF65412E65401843C95729C4CFA2FA3A773C8AF928884F5AC37F0A088D000758301A811009C68F30B72F4F6A15EC9FA429E81C8E912E308434AD064658B158076E4B7AF526EBE25B78E6549E3201FF2F00715230033A1A6B54FA2D4001D7F102D90838139A91C05A6A9D68F882B7DE89244F65592F05B20F3FB377B7E44DC44891392192A62BF1B9D54DC548072BFA6F9C0BE264C1352A0D1737834831BDAD11526E28B98BC6C736BBAB1468F972A738BA451261E2151E17C694370D82615B7ED5A6147F0182B520E6B15768F7AF395453B5BE1CB67E7EB14D346DA7AC90108001A10071C3FFA75C404737C2D1C0698DDDBCCD06C4F596950EA6D31EB87CBC914E1BACFA7590A02D529F902BE4826DD4530ECCAD7D346CABA951E897718027047297FA7A4065CD6631FEB621BDB69F3FE937EC4B7A46A957A055F6EEDFB8A4C1E4ABBF0A8900A5636EF45A2DC60BDB448613282EF750BAFB67ED87F4F88BD717C821C20FE91509E128A034915AEAF8AEC9D67C171EEF164A9062B155C1D445003C3140F10E177F1B0671079B91D045AEC8D89A312567D9857D5F6F01E1CC9BCBDF9EDC38C2A2015EC4E0502713243B5AE5EA87D286FCF9CB08BEA29DB65F64300CBD2BC7528A5D513D43ECA03DC9E6C07781C529BA2042E1476682A194F5D83884CD5763CAD47F67C822772DD435AB812D7F1F0F11F9A592D5497144A7AE497AC53736AA4E5F6C3B40850D51688625E53669F8A9F9B14813FB99D1C6DE164C2767D181C06FC229CEA6EAC4AF180176DAD5132FB2BFCA6AFFBA4E08F39F97509E1D6B2C4EAFC98FED717DB42688A117D6FDE5E72342E65C787B8FF305E5C5DA775CAFBC66A5CEB55C3A72CC516CFD6B6BE44222EC251F76B22E0F6336416FC412BF2F28600FC6CFA2FEBE13A6D36DD99F247BDF1B5F90BED48CF458FF864628BEF87C1753B6BD63C058DE888702EA5E732E508CD3A39318214E1499F2D82DB19951E19361A8015603F2D9EC7E0EC0139AB910603D255728B708484D434CDEBFCB5CE8F63C312BEC8E9370F28972B9E30096E17443C23FF92412D35CD7797A5B2439D6C78200F4EC5B673553F4EECE004F26D0836DFD0BE02312DF38BA25AA06B8AE4BEF06086CD29F34C98C79831C778DBE6BDF73A5EBAFAD1C71986E81F7D0D70AB970870ABA617E186FB7165FAC1550E7D29F7BB1264E8C02A0A99F199E3CD90BD766CC5FEBC7D054754220EF18D8AEB8DC745FE1AC8038F0EC6D35A34EB89D1422727272085E712CD29E9F60CE00D88D35214A439C34F79CAFEBE4B8B9A9191C6C89550C982CBEFD98208207AD46CDCAF23FF3EFA244B3699BA0DADD608543E8BD8AA6FCEC3CF0819DA1269FE72F70604D941EF4336889EC98E29AE03F02D2173EC7DA1325E65E89EFDBC4C5463281413A958FAE77860EEFA761B57B8FC07FDA0807C60A91834228B0FF85A314CF0FAE6B7406354D689489C6037953DBD19AD0A8B56C125B02D0D88AA02577C1F977DCBE01DEAA2863F108A1DF24DF3C6C5C6256E35495B41B29722AD388B6E09C25E9161CACBEAE1DC260B92FE8C2D1048F9AA59D8C60E2CDF72B9D8DD6AD09237D04DBFAB40640303838BBA8F2446703597C7E10072557603849CE7D030911C99213FF5FEB9C6158C563A640D754F62486ECF3516A0D406343C558BBE03E6F721AC231315029901A072EF87B266D5C49903A15E8C7DC9BA654F38A15346B1CE9C3968829AF2D925CE4DDC0A0BFE540BD1D2037FC422A157D1D34F8B3CE083E48826EDF09822606D2AF227D7B602249A3915B5A2A03240B6942CB535AAF2938B0CAE47C2888B52D6E56891C029A73176B11E80B7B8D692C970DFC40F74C389B459E20C486FAF03ABADA84D80225E3D8B5FD37E3BD7B9F52267BD0EDDF370370B";
        let secret_key = SecretKey::decode(&hex_to_bytes(secret_key_binary).unwrap());
        let trapdoor_output =
            TrapdoorOutput::decode(&hex_to_bytes(trapdoor_output_binary).unwrap());

        let festa = FESTACrypto::new();
        let inverted = festa.trapdoor_inverse(&trapdoor_output, &secret_key);
        println!("inverted : {}", inverted);
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
