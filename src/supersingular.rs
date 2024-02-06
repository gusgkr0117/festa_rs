//! Functions for the supersingular curves
use crate::{ecFESTA::{factored_kummer_isogeny, Curve, KummerLineIsogeny, Point}, fields::FpFESTAExt::{BASIS_ORDER, L_POWER}, pairing::tate_pairing};
use num_bigint::BigUint;


/// In the case where we remove a known large even factor, it's faster to do this first
/// to save on addition in xDBLADD
pub fn clear_cofactor(E: &Curve, P: &Point, k : &BigUint, even_power: usize) -> Point {
    let Q = E.double_iter(&P, even_power);
    let x = k >> even_power;
    E.mul_big(&Q, &x)
}

/// Generate a random D-torsion point
pub fn generate_point_order_factored_D(E: &Curve, D: &BigUint, even_power: usize) -> Point {
    let mut rng = rand::thread_rng();
    let cofactor = BigUint::from_slice(&BASIS_ORDER) / D;

    loop {
        let mut random_point = E.rand_point(&mut rng);
        random_point = clear_cofactor(E, &random_point, &cofactor, even_power);
        
        if random_point.isinfinity() != 0 {
            continue;
        }

        if E.mul_big(&random_point, D).isinfinity() == 0 {
            continue;
        }

        return random_point;
    }
}

/// Generate a basis of E(Fq)[D] of supersingular curve
pub fn torsion_basis(E: &Curve, D: &BigUint, even_power: usize) -> (Point, Point) {
    let P = generate_point_order_factored_D(E, D, even_power);
    loop {
        let Q = generate_point_order_factored_D(E, D, even_power);
        let ePQ = tate_pairing(E, &P, &Q, D);
        // TODO : check if the order of ePQ is equal to D
        return (P, Q)
    }
}

/// Generate a random isogeny
pub fn random_isogeny_x_only(E: &Curve, D: &[(u32, u32)], steps: usize) -> (Vec<KummerLineIsogeny>, Curve) {
    let mut isogeny_list = Vec::new();
    let one_step_degree = D.to_owned();
    let mut degree = BigUint::from(1u32);
    let mut domain_curve = E.clone();

    for (l, e) in one_step_degree.iter() {
        degree = degree * BigUint::from(*l).pow(*e);
    }

    for _ in 0..steps {
        let kernel_point = generate_point_order_factored_D(&domain_curve, &degree, L_POWER as usize);
        let mut phi = factored_kummer_isogeny(&domain_curve, &kernel_point, &one_step_degree);
        domain_curve = phi.last().expect("phi is empty").get_codomain();
        isogeny_list.append(&mut phi);
    }

    (isogeny_list, domain_curve)
}

#[cfg(test)]
mod tests {
    use crate::ecFESTA::Fq;

    use super::*;
    #[test]
    fn compute_random_isogeny() {
        let start_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let degree : [(u32, u32); 3] = [(59, 1), (3359, 1), (6299, 1)];
        let (phi, coE) = random_isogeny_x_only(&start_curve, &degree, 1);
        println!("codomain curve : {}", coE.j_invariant());
        for isog in phi.iter() {
            println!("{:#}", isog);
        }
    }
}