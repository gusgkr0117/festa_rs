/// Implementation of solving DLP in the finite field Fq
/// Use Baby-Step Giant-Step algorithm and Pohlig-Hellman algorithm
use crate::{
    ecFESTA::{Curve, Fq, Point},
    pairing::weil_pairing,
};
use num_bigint::{BigInt, BigUint, Sign, ToBigInt};

/// Computing the extended GCD algorithm of given two BigUint type numbers
/// compute r such that a*r - b*y = d = gcd(a,b)
pub fn xgcd_big(a: &BigUint, b: &BigUint) -> BigUint {
    let (mut s1, mut t1) = (BigInt::from(1u32), BigInt::from(0u32));
    let (mut s2, mut t2) = (BigInt::from(0u32), BigInt::from(1u32));
    let (mut m, mut n) = (a.to_bigint().unwrap(), b.to_bigint().unwrap());

    if m < n {
        (m, n) = (n, m);
        ((s1, t1), (s2, t2)) = ((s2, t2), (s1, t1));
    }

    while n != BigInt::from(0) {
        let q = &m / &n;
        let r = &m - &q * &n;
        (m, n) = (n, r);
        ((s1, t1), (s2, t2)) = ((s2.clone(), t2.clone()), (s1 - &q * s2, t1 - &q * t2));
    }

    s1 = s1 % b.to_bigint().unwrap();

    if s1.sign() == Sign::Minus {
        s1 += b.to_bigint().unwrap();
    }

    s1.to_biguint().unwrap()
}

/// Computing Chinese Remainder Theorem(CRT) for given BigUint type numbers as inputs
pub fn compute_crt(remainder_list: &Vec<BigUint>, modulus_factor: &[(u32, u32)]) -> BigUint {
    let mut result = BigUint::from(0u32);
    let modulus = modulus_factor
        .iter()
        .fold(BigUint::from(1u32), |s, (l, e)| {
            s * BigUint::from(*l).pow(*e)
        });

    for (i, (l, e)) in modulus_factor.iter().enumerate() {
        let factor = BigUint::from(*l).pow(*e);
        let mut cofactor = BigUint::from(1u32);
        for (l1, e1) in modulus_factor.iter() {
            if *l1 != *l {
                cofactor *= BigUint::from(*l1).pow(*e1);
            }
        }
        let cofactor_inv = xgcd_big(&cofactor, &factor);
        result = (result + &remainder_list[i] * cofactor_inv * cofactor) % &modulus;
    }
    result
}

/// Computing the square root of a small integer of u32 type
pub fn sqrt_floor(n: u32) -> u32 {
    if n == 1 {
        return 1u32;
    }

    let (mut f, mut b) = (0u32, n);
    while f + 1 < b {
        let x = (f + b) / 2;
        let sqr = (x as u64) * (x as u64);
        if sqr > (n as u64) {
            b = x;
        } else if sqr < (n as u64) {
            f = x;
        } else {
            return x;
        }
    }

    if b * b == n {
        return b;
    }

    f
}

/// Baby-Step Giant-Step(BSGS) algorithm for a small prime order
pub fn bsgs(beta: &Fq, alpha: &Fq, order: u32) -> u32 {
    let m = sqrt_floor(order) + 1;
    let mut alpha_list = Vec::new();
    let mut gamma = Fq::ONE;
    for i in 0..m {
        alpha_list.push((i, gamma));
        gamma *= alpha;
    }

    let alpha_m_inv = gamma.invert();
    let mut target = *beta;
    for i in 0..m {
        if let Some((j, _)) = alpha_list.iter().find(|(_, v)| v.equals(&target) != 0) {
            return (i * m + j) % order;
        }

        target *= alpha_m_inv;
    }
    panic!("No bsgs solution");
}

/// Solving a DLP problem using Pohlig-Hellman algorithm
/// where the order of the group is a power of a prime : l^e
pub fn ph_dlp_power(beta: &Fq, alpha: &Fq, l: u32, e: u32) -> BigUint {
    let mut h;
    let mut x = BigUint::from(0u32);
    let mut gamma = *alpha;
    let mut l_i = BigUint::from(1u32);
    for _ in 0..e - 1 {
        gamma = gamma.pow_small(l);
    }

    debug_assert!(
        gamma.equals(&Fq::ONE) == 0 && gamma.pow_small(l).equals(&Fq::ONE) != 0,
        "The order of gamma is wrong"
    );

    for i in 0..e {
        let g_x_inv = alpha.pow_big(&x).invert();
        h = *beta * g_x_inv;
        for _ in 0..e - 1 - i {
            h = h.pow_small(l);
        }
        debug_assert!(
            h.pow_small(l).equals(&Fq::ONE) != 0,
            "The order of h is wrong"
        );

        let d = bsgs(&h, &gamma, l);
        debug_assert!(gamma.pow_small(d).equals(&h) != 0, "bsgs is wrong");
        x = x + &l_i * d;
        l_i *= l;
    }

    x
}

/// Solving a DLP problem using Pohlig-Hellman algorithm
/// Get a factored order as an input to efficiently compute the solution
/// Compute x such that beta = alpha^x in Fq where alpha is an element of the given order
pub fn ph_dlp(beta: &Fq, alpha: &Fq, factored_order: &[(u32, u32)]) -> BigUint {
    let mut x = Vec::new();
    for (l, e) in factored_order.iter() {
        let mut alpha_i = *alpha;
        let mut beta_i = *beta;
        for (l1, e1) in factored_order.iter() {
            if *l1 != *l {
                for _ in 0..*e1 {
                    alpha_i = alpha_i.pow_small(*l1);
                    beta_i = beta_i.pow_small(*l1);
                }
            }
        }
        let x_i = ph_dlp_power(&beta_i, &alpha_i, *l, *e);
        x.push(x_i);
    }

    compute_crt(&x, factored_order)
}

/// Given a basis P, Q of E[D] finds
/// a, b such that R = [a]P + [b]Q.
pub fn bidlp(
    E: &Curve,
    R: &Point,
    P: &Point,
    Q: &Point,
    factored_D: &[(u32, u32)],
    ePQ: Option<Fq>,
) -> (BigUint, BigUint) {
    let order = factored_D.iter().fold(BigUint::from(1u32), |r, (l, e)| {
        r * BigUint::from(*l).pow(*e)
    });
    let pair_PQ = match ePQ {
        Some(x) => x,
        None => weil_pairing(E, P, Q, &order),
    };
    let pair_a = weil_pairing(E, R, Q, &order);
    let pair_b = weil_pairing(E, R, &-P, &order);
    let a = ph_dlp(&pair_a, &pair_PQ, factored_D);
    let b = ph_dlp(&pair_b, &pair_PQ, factored_D);

    (a, b)
}

#[cfg(test)]
mod tests {
    use crate::{
        ecFESTA::Curve,
        fields::FpFESTAExt::{D1, D1_FACTORED, L_POWER},
        pairing::tate_pairing,
        supersingular::torsion_basis,
    };

    use super::*;
    #[test]
    fn compute_dlp_using_pairing() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let (P, Q) = torsion_basis(&test_curve, &D1_FACTORED, L_POWER as usize);
        let d1 = BigUint::from_slice(&D1);
        let ePQ = tate_pairing(&test_curve, &P, &Q, &d1);
        let ePnQ = tate_pairing(&test_curve, &P, &test_curve.mul_small(&Q, 19u64), &d1);

        println!("ePQ : {}", ePQ.pow_small(19u32));
        println!("ePnQ : {}", ePnQ);

        let dlp = ph_dlp(&ePnQ, &ePQ, &D1_FACTORED);
        println!("dlp result : {}", dlp);
    }

    #[test]
    fn compute_xgcd() {
        let n = BigUint::from(1241u32);
        let m = BigUint::from(512u32);
        let r = xgcd_big(&m, &n);
        println!("r : {}", r);
        println!("(m * r) % n : {}", (m * r) % n);
    }

    #[test]
    fn test_crt() {
        let modulus_list: [(u32, u32); 3] = [(3, 2), (5, 1), (7, 3)];
        let remainder_list: Vec<BigUint> = vec![
            BigUint::from(5u32),
            BigUint::from(3u32),
            BigUint::from(23u32),
        ];
        let result = compute_crt(&remainder_list, &modulus_list);
        println!("result : {}", result);
        println!("result % 3^2 = {}", &result % BigUint::from(9u32));
        println!("result % 5^1 = {}", &result % BigUint::from(5u32));
        println!("result % 7^3 = {}", &result % BigUint::from(7u32).pow(3));
    }

    #[test]
    fn test_bidlp() {
        let factored_D = [(71u32, 3u32)];
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let (P, Q) = torsion_basis(&test_curve, &factored_D, L_POWER as usize);
        let R = test_curve.add(
            &test_curve.mul_small(&P, 132u64),
            &test_curve.mul_small(&Q, 3142u64),
        );
        let (a, b) = bidlp(&test_curve, &R, &P, &Q, &factored_D, None);

        let test_point = test_curve.add(&test_curve.mul_big(&P, &a), &test_curve.mul_big(&Q, &b));

        println!("R = [{}]*P + [{}]*Q", a, b);
        assert!(R.equals(&test_point) != 0, "bidlp is wrong");
    }
}
