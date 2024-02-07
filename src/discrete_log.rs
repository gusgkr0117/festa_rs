/// Implementation of solving DLP in the finite field Fq
/// Use Baby-Step Giant-Step algorithm and Pohlig-Hellman algorithm
use crate::ecFESTA::Fq;
use num_bigint::BigUint;

/// Computing the extended GCD algorithm of given two BigUint type numbers
/// compute r such that a*r - b*y = d = gcd(a,b)
pub fn xgcd_big(a: &BigUint, b: &BigUint) -> BigUint {
    todo!("implementation of extended GCD algorithm for BigUints");
    BigUint::from(1u32)
}

/// Computing Chinese Remainder Theorem(CRT) for given BigUint type numbers as inputs
pub fn compute_crt(remainder_list: &Vec<BigUint>, modulus_factor: &[(u32, u32)]) -> BigUint {
    let mut result = BigUint::from(0u32);

    for (i, (l, e)) in modulus_factor.iter().enumerate() {
        let factor = BigUint::from(*l).pow(*e);
        let mut cofactor = BigUint::from(1u32);
        for (l1, e1) in modulus_factor.iter() {
            if *l1 != *l {
                cofactor *= BigUint::from(*l1).pow(*e1);
            }
        }
        let cofactor_inv = xgcd_big(&cofactor, &factor);
        result += &remainder_list[i] * cofactor_inv * cofactor;
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

/// Baby-step and Giant-step algorithm for a small prime order
pub fn bsgs(beta: &Fq, alpha: &Fq, order: u32) -> u32 {
    let m = sqrt_floor(order);
    let mut alpha_list = Vec::new();
    let mut gamma = Fq::ONE;
    for i in 0..m {
        alpha_list.push((i, gamma));
        gamma *= alpha;
    }

    let alpha_m_inv = (gamma * alpha).invert();
    let mut target = *beta;
    for i in 0..m {
        if let Some((j, _)) = alpha_list.iter().find(|(_, v)| v.equals(&target) != 0) {
            return i * m + j;
        }

        target *= alpha_m_inv;
    }
    panic!("No bsgs solution");
}

/// Solving a DLP problem using Pohlig-Hellman algorithm
/// where the order of the group is a power of a prime : l^e
pub fn ph_dlp_power(beta: &Fq, alpha: &Fq, l: u32, e: u32) -> BigUint {
    let mut h = beta.pow_small(l);
    let mut x = BigUint::from(0u32);
    let mut gamma = *alpha;
    let mut l_i = BigUint::from(1u32);
    for _ in 0..e - 1 {
        gamma = gamma.pow_small(l);
    }

    for i in 0..e {
        h = *beta;
        for _ in 0..e - 1 - i {
            h = h.pow_small(l);
        }
        let d = bsgs(&h, &gamma, l);
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn compute_dlp_using_pairing() {}
}
