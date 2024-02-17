//! Functions for the supersingular curves
use std::collections::HashMap;

use crate::{
    ecFESTA::{factored_kummer_isogeny, Curve, Fq, KummerLineIsogeny, Point},
    fields::{
        FpFESTA::Fp,
        FpFESTAExt::{BASIS_ORDER, BIGUINT_MODULUS, L_POWER},
    },
    pairing::weil_pairing,
};
use num_bigint::BigUint;

/// Given the factored order D as an input, check if the point has an order D
pub fn point_has_factored_order(E: &Curve, P: &Point, factored_order: &[(u32, u32)]) -> bool {
    let mut D_top = BigUint::from(1u32);
    let mut pis = Vec::new();

    for (l, e) in factored_order.iter() {
        D_top *= BigUint::from(*l).pow(e - 1);
        pis.push(BigUint::from(*l));
    }

    let G_top = E.mul_big(P, &D_top);
    if G_top.equals(&Point::INFINITY) != 0 {
        return false;
    }

    /// Input : A list of elements 'G_list', such that G is the first entry and the rest is empty
    /// in the sublist G_list[lower:upper]
    fn batch_cofactor_mul(
        E: &Curve,
        G_list: &mut Vec<Point>,
        pis: &Vec<BigUint>,
        lower: usize,
        upper: usize,
    ) {
        if lower > upper {
            panic!("Wrong input to cofactor_multiples()");
        }

        if upper - lower == 1 {
            return;
        }

        let mid = lower + (upper - lower + 1) / 2;
        let (mut cl, mut cu) = (BigUint::from(1u32), BigUint::from(1u32));
        for i in lower..mid {
            cu *= &pis[i as usize];
        }

        for i in mid..upper {
            cl *= &pis[i as usize];
        }

        G_list[mid] = E.mul_big(&G_list[lower], &cu);
        G_list[lower] = E.mul_big(&G_list[lower], &cl);

        batch_cofactor_mul(E, G_list, pis, lower, mid);
        batch_cofactor_mul(E, G_list, pis, mid, upper);
    }

    let mut G_list = Vec::new();
    for _ in 0..pis.len() {
        G_list.push(Point::INFINITY);
    }

    G_list[0] = G_top;
    if pis.len() > 1 {
        batch_cofactor_mul(E, &mut G_list, &pis, 0, pis.len());
        return G_list.iter().fold(true, |s, x| s & (x.isinfinity() == 0));
    }

    true
}

/// Given the factored order D as an input, check if the group element has an order D
pub fn has_factored_order(v: Fq, factored_order: &[(u32, u32)]) -> bool {
    let mut D_top = BigUint::from(1u32);
    let mut pis = Vec::new();

    for (l, e) in factored_order.iter() {
        D_top *= BigUint::from(*l).pow(e - 1);
        pis.push(BigUint::from(*l));
    }

    let G_top = v.pow_big(&D_top);
    if G_top.equals(&Fq::ONE) != 0 {
        return false;
    }

    /// Input : A list of elements 'G_list', such that G is the first entry and the rest is empty
    /// in the sublist G_list[lower:upper]
    fn batch_cofactor_mul(G_list: &mut Vec<Fq>, pis: &Vec<BigUint>, lower: usize, upper: usize) {
        if lower > upper {
            panic!("Wrong input to cofactor_multiples()");
        }

        if upper - lower == 1 {
            return;
        }

        let mid = lower + (upper - lower + 1) / 2;
        let (mut cl, mut cu) = (BigUint::from(1u32), BigUint::from(1u32));
        for i in lower..mid {
            cu *= &pis[i as usize];
        }

        for i in mid..upper {
            cl *= &pis[i as usize];
        }

        G_list[mid] = G_list[lower].pow_big(&cu);
        G_list[lower] = G_list[lower].pow_big(&cl);

        batch_cofactor_mul(G_list, pis, lower, mid);
        batch_cofactor_mul(G_list, pis, mid, upper);
    }

    let mut G_list = Vec::new();
    for _ in 0..pis.len() {
        G_list.push(Fq::ONE);
    }

    G_list[0] = G_top;
    if pis.len() > 1 {
        batch_cofactor_mul(&mut G_list, &pis, 0, pis.len());

        return G_list
            .iter()
            .fold(true, |s, x| s & (x.equals(&Fq::ONE) == 0));
    }

    true
}

/// In the case where we remove a known large even factor, it's faster to do this first
/// to save on addition in xDBLADD
pub fn clear_cofactor(E: &Curve, P: &Point, k: &BigUint, even_power: usize) -> Point {
    let Q = E.double_iter(&P, even_power);
    let x = k >> even_power;
    E.mul_big(&Q, &x)
}

/// Generate a random D-torsion point canonically
pub fn generate_point_order_factored_D(
    E: &Curve,
    factored_order: &[(u32, u32)],
    even_power: usize,
) -> Point {
    let order = factored_order
        .iter()
        .fold(BigUint::from(1u32), |m, (l, e)| {
            m * BigUint::from(*l).pow(*e)
        });
    let cofactor = BigUint::from_slice(&BASIS_ORDER) / &order;

    for random_point in E.iter() {
        let torsion_point = clear_cofactor(E, &random_point, &cofactor, even_power);

        if torsion_point.isinfinity() != 0 {
            continue;
        }

        if !point_has_factored_order(E, &torsion_point, factored_order) {
            continue;
        }

        return torsion_point;
    }

    Point::INFINITY
}

/// Generate a basis of E(Fq)[D] of supersingular curve
pub fn torsion_basis(E: &Curve, factored_D: &[(u32, u32)], even_power: usize) -> (Point, Point) {
    let D = factored_D.iter().fold(BigUint::from(1u32), |m, (l, e)| {
        m * BigUint::from(*l).pow(*e)
    });
    let cofactor = BigUint::from_slice(&BASIS_ORDER) / &D;

    let P = generate_point_order_factored_D(E, factored_D, even_power);
    for random_point in E.iter() {
        let Q = clear_cofactor(E, &random_point, &cofactor, even_power);
        if P.equals(&Q) != 0 || Q.isinfinity() != 0 {
            continue;
        }
        if E.mul_big(&Q, &D).isinfinity() == 0 {
            continue;
        }

        if !point_has_factored_order(E, &Q, factored_D) {
            continue;
        }

        let ePQ = weil_pairing(E, &P, &Q, &D);
        //check if the order of ePQ is equal to D
        if !has_factored_order(ePQ, factored_D) {
            continue;
        }
        return (P, Q);
    }

    (Point::INFINITY, Point::INFINITY)
}

/// Optimized algorithm for generating k*2^l torsion basis
pub fn entangled_torsion_basis(E: &Curve, cofactor: &BigUint) -> (Point, Point) {
    fn precompute_elligator_tables() -> (Vec<(Fq, Fq)>, Vec<(Fq, Fq)>) {
        let u = Fq::ZETA * Fq::TWO;
        let mut T1 = Vec::<(Fq, Fq)>::new();
        let mut T2 = Vec::<(Fq, Fq)>::new();
        let mut r = Fq::ONE;
        for _ in 0..30 {
            let v = (Fq::ONE + &u * &r.square()).invert();
            if v.legendre() == 1 {
                T2.push((r, v));
            } else {
                T1.push((r, v));
            }
            r += Fq::ONE;
        }

        (T1, T2)
    }

    let p = BigUint::from_slice(&BIGUINT_MODULUS);
    let p_sqrt = (p + BigUint::from(1u32)) / BigUint::from(4u32);
    let u = Fq::ZETA * Fq::TWO;
    let u0 = Fq::ONE + Fq::ZETA;

    let (TQNR, TQR) = precompute_elligator_tables();

    // Pick the lookup table depending on whether
    // A = a + ib is a quadratic residue or quadratic nonresidue
    let A = E.get_constant();
    let T = match A.legendre() {
        1i32 => TQNR,
        -1i32 => TQR,
        _ => panic!("invalid curve for entangled torsion basis algorithm"),
    };

    // Look through the table to find point with
    // rational (x, y)
    let mut y: Option<Fq> = None;
    let (mut s, mut c, mut d): (Fp, Fp, Fp) = (Fp::ONE, Fp::ONE, Fp::ONE);
    let (mut x, mut r) = (Fq::ONE, Fq::ONE);
    for (ri, v) in T.iter() {
        (r, x) = (*ri, -A * v);
        let t = x * (x.square() + A * x + Fq::ONE);

        // break when we find rational y: t = y^2
        (c, d) = t.get_coeff();
        let z = c.square() + d.square();
        s = z.pow_big(&p_sqrt);
        if s.square().equals(&z) != 0 {
            y = Some(t.sqrt().0);
            break;
        }
    }

    if y.is_none() {
        panic!("Never found a y-coordinate, increase the lookup table size");
    }

    let z = (&c + &s).half();
    let alpha = z.pow_big(&p_sqrt);
    let beta = d / (&alpha + &alpha);

    let y = match alpha.square().equals(&z) != 0 {
        true => Fq::new(&alpha, &beta),
        false => -Fq::new(&beta, &alpha),
    };

    let S1 = Point::new_xy(&x, &y);

    let (s21, s22) = ((u * r.square() * x), (u0 * r * y));

    let S2 = Point::new_xy(&s21, &s22);

    (E.mul_big(&S1, &cofactor), E.mul_big(&S2, &cofactor))
}

/// Generate a random isogeny
pub fn random_isogeny_x_only(
    E: &Curve,
    D: &[(u32, u32)],
    steps: usize,
) -> (Vec<KummerLineIsogeny>, Curve) {
    let mut isogeny_list = Vec::new();
    let one_step_degree = D.to_owned();
    let mut degree = BigUint::from(1u32);
    let mut domain_curve = E.clone();

    for (l, e) in one_step_degree.iter() {
        degree = degree * BigUint::from(*l).pow(*e);
    }

    for _ in 0..steps {
        let kernel_point = generate_point_order_factored_D(&domain_curve, &D, L_POWER as usize);
        let mut phi = factored_kummer_isogeny(&domain_curve, &kernel_point, &one_step_degree);
        domain_curve = phi.last().expect("phi is empty").get_codomain();
        isogeny_list.append(&mut phi);
    }

    (isogeny_list, domain_curve)
}

/// Computes a D "odd" degree isogeny from E using
/// x-only arithmetic and returns the KummerIsogeny
/// together with the codomain curve
///
/// The kernel of this isogeny is K = <P + [m]Q>
/// where P, Q are the canonical D torsion basis
pub fn isogeny_from_scalar_x_only(
    E: &Curve,
    factored_D: &[(u32, u32)],
    m: &BigUint,
    basis: Option<(Point, Point)>,
) -> (Vec<KummerLineIsogeny>, Curve) {
    let (P, Q) = match basis {
        Some(x) => x,
        None => torsion_basis(E, factored_D, L_POWER as usize),
    };

    let K = E.add(&P, &E.mul_big(&Q, &m));
    let phi = factored_kummer_isogeny(&E, &K, factored_D);
    let codomain = phi.last().unwrap().get_codomain();

    (phi, codomain)
}

#[cfg(test)]
mod tests {
    use crate::{
        ecFESTA::Fq,
        fields::FpFESTAExt::{D1, D1_FACTORED, D2, D2_FACTORED},
        pairing::weil_pairing,
    };

    use super::*;
    #[test]
    fn compute_random_isogeny() {
        let start_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let degree: [(u32, u32); 3] = [(59, 1), (3359, 1), (6299, 1)];
        let (phi, coE) = random_isogeny_x_only(&start_curve, &degree, 1);
        println!("codomain curve : {}", coE.j_invariant());
        for isog in phi.iter() {
            println!("{:#}", isog);
        }
    }

    fn compute_torsion_point(d: &BigUint, factored_d: &[(u32, u32)]) {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let P = generate_point_order_factored_D(&test_curve, &factored_d, L_POWER as usize);
        println!("d torsion point : {}", P);
        assert!(
            test_curve.mul_big(&P, &d).isinfinity() != 0,
            "[d] * P is nonzero"
        );
        for (l, _) in factored_d.iter() {
            let cofactor = d / &BigUint::from(*l);
            assert!(
                test_curve.mul_big(&P, &cofactor).isinfinity() == 0,
                "[d/{}] * P is zero",
                l
            );
        }
    }

    #[test]
    fn compute_torsion_d1_and_d2() {
        compute_torsion_point(&BigUint::from_slice(&D1), &D1_FACTORED);
        compute_torsion_point(&BigUint::from_slice(&D2), &D2_FACTORED);
    }

    #[test]
    fn compute_torsion_basis() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));

        let (P, Q) = torsion_basis(&test_curve, &[(2, L_POWER)], 0 as usize);

        assert!(point_has_factored_order(&test_curve, &P, &[(2, L_POWER)]));
        assert!(point_has_factored_order(&test_curve, &Q, &[(2, L_POWER)]));

        let l_power = BigUint::from(2u32).pow(L_POWER);
        let ePQ = weil_pairing(&test_curve, &P, &Q, &l_power);

        println!("e(P,Q)^r = {}", ePQ.pow_big(&l_power).pow_small(1024u32));

        // println!("P = {P}");
        // println!("Q = {Q}");
    }

    #[test]
    fn compute_l_power_torsion_basis() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let basis_order = BigUint::from_slice(&BASIS_ORDER);
        let cofactor = &basis_order / &l_power;
        let (P, Q) = entangled_torsion_basis(&test_curve, &cofactor);
        println!("P : {}", P);
        println!("Q : {}", Q);

        assert_ne!(
            test_curve.mul_big(&P, &l_power).isinfinity(),
            0,
            "P has wrong order"
        );
        assert_ne!(
            test_curve.mul_big(&Q, &basis_order).isinfinity(),
            0,
            "Q has wrong order"
        );

        let cofactor = BigUint::from(2u32).pow(L_POWER - 11);
        let (R, S) = (
            test_curve.mul_big(&P, &cofactor),
            test_curve.mul_big(&Q, &cofactor),
        );
        // Test if two given points are independent
        for (i, j) in (1..2048).zip(1..2048) {
            if test_curve
                .mul_small(&R, i)
                .equals(&test_curve.mul_small(&S, j))
                != 0
            {
                println!("{j}S = {i}R");
                return;
            }
        }
        println!("R and S are linearly independent");
    }

    #[test]
    fn iterate_curve() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        for (i, point) in test_curve.iter().enumerate() {
            println!("point : {}", point);
            if i > 5 {
                break;
            }
        }
    }
}
