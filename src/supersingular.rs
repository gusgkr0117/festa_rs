//! Functions for the supersingular curves
use crate::{
    ecFESTA::{factored_kummer_isogeny, Curve, Fq, KummerLineIsogeny, Point},
    fields::FpFESTAExt::{BASIS_ORDER, L_POWER},
    pairing::tate_pairing,
};
use num_bigint::BigUint;

/// Given the factored order D as an input, check if the group element has D
pub fn has_factored_order(v: Fq, factored_order: &[(u32, u32)]) -> bool {
    let mut D_top = BigUint::from(1u32);
    let mut pis = Vec::new();

    for (l, e) in factored_order.iter() {
        D_top *= BigUint::from(*l).pow(e - 1);
        pis.push(BigUint::from(*l));
    }

    let G_top = v.pow_big(&D_top);
    if G_top.iszero() != 0 {
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
        let mut sum = true;
        for element in G_list.iter().map(|x| x.iszero() != 0).into_iter() {
            sum &= !element;
        }
        return sum;
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
pub fn torsion_basis(E: &Curve, factored_D: &[(u32, u32)], even_power: usize) -> (Point, Point) {
    let mut D = BigUint::from(1u32);
    for (l, e) in factored_D.iter() {
        D *= BigUint::from(*l).pow(*e);
    }

    let P = generate_point_order_factored_D(E, &D, even_power);
    loop {
        let Q = generate_point_order_factored_D(E, &D, even_power);
        let ePQ = tate_pairing(E, &P, &Q, &D);
        //check if the order of ePQ is equal to D
        if !has_factored_order(ePQ, factored_D) {
            continue;
        }
        return (P, Q);
    }
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
        let kernel_point =
            generate_point_order_factored_D(&domain_curve, &degree, L_POWER as usize);
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
        let degree: [(u32, u32); 3] = [(59, 1), (3359, 1), (6299, 1)];
        let (phi, coE) = random_isogeny_x_only(&start_curve, &degree, 1);
        println!("codomain curve : {}", coE.j_invariant());
        for isog in phi.iter() {
            println!("{:#}", isog);
        }
    }

    #[test]
    fn compute_torsion_basis() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let degree: [(u32, u32); 3] = [(59, 1), (3359, 1), (6299, 1)];
        let (P, Q) = torsion_basis(&test_curve, &degree, L_POWER as usize);
        println!("P : {}", P);
        println!("Q : {}", Q);
    }
}
