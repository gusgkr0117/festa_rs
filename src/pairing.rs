use std::ops::Mul;

use num_bigint::BigUint;

use crate::thetaFESTA::{Curve, Fq, Point};

/// Evalute the line g crossing the two points P and Q
/// and output g(R)
fn eval_line(E: &Curve, P: &Point, Q: &Point, R: &Point) -> Fq {
    let result;

    let (Px, Py) = P.to_xy();
    let (Qx, Qy) = Q.to_xy();
    let (Rx, Ry, Rz) = R.to_xyz();
    let a = E.get_constant();
    // When P or Q is the point at infinity i.e. D(g) = 0,
    // the dividing by this evaluation of the line must be ignored
    debug_assert!(P.isinfinity() == 0 && Q.isinfinity() == 0);

    if P.equals(&Q) != 0 {
        if Py.iszero() != 0 {
            result = Rx + (-Px) * Rz;
            return result;
        }

        let lambda = (Px.square().mul3() + (a * Px).mul2() + Fq::ONE) / Py.mul2();
        result = -lambda * Rx + Ry + (-Py + lambda * Px) * Rz;
    } else if P.equals(&(-Q)) != 0 {
        result = Rx + (-Px) * Rz;
    } else {
        let lambda = (Qy - Py) / (Qx - Px);
        result = -lambda * Rx + Ry + (-Py + lambda * Px) * Rz;
    }
    result
}

/// Evaluate Tate pairing using Miller's algorithm
/// If two points P and Q are linearly dependent, it outputs zero.
/// Reference : [crypto.stanford.edu](https://crypto.stanford.edu/pbc/notes/ep/miller.html)
pub fn tate_pairing(E: &Curve, P: &Point, Q: &Point, order: &BigUint) -> Fq {
    let mut f = Fq::ONE;
    let mut V = P.clone();

    for t in (0..order.bits() - 1).rev() {
        let S = E.double(&V);
        if S.isinfinity() == 0 {
            f = f.square() * eval_line(E, &V, &V, Q) / eval_line(E, &S, &-S, Q);
        } else {
            let (_, _, Qz) = Q.to_xyz();
            f = f.square() * eval_line(E, &V, &V, Q) / Qz;
        }
        V = S;
        if order.bit(t) {
            let S = E.add(&V, P);
            if S.isinfinity() == 0 {
                f = f * eval_line(E, &V, P, Q) / eval_line(E, &S, &-S, Q);
            } else {
                let (_, _, Qz) = Q.to_xyz();
                f = f * eval_line(E, &V, P, Q) / Qz;
            }
            V = S;
        }
    }

    debug_assert!(V.isinfinity() != 0, "pairing point is invalid");

    if f.iszero() != 0 {
        return Fq::ONE;
    }

    f
}

/// Computing Weil pairing of given two points using Tate pairing(Miller's algorithm).
pub fn weil_pairing(E: &Curve, P: &Point, Q: &Point, order: &BigUint) -> Fq {
    let ePQ = tate_pairing(E, P, Q, order);
    let eQP = tate_pairing(E, Q, P, order);

    if order.bit(0) {
        -ePQ / eQP
    } else {
        ePQ / eQP
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use crate::{
        fields::FpFESTAExt::{BASIS_ORDER, D1, D1_FACTORED, L_POWER},
        supersingular::{entangled_torsion_basis, has_factored_order, torsion_basis},
    };

    use super::*;
    use num_bigint::BigUint;
    #[test]
    fn compute_pairing() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let basis_order = BigUint::from_slice(&BASIS_ORDER);

        // sample order for test
        let new_order = BigUint::from_slice(&D1);

        assert!((&basis_order % &new_order) == BigUint::from(0u32));

        let (test_P, test_Q) = torsion_basis(&test_curve, &D1_FACTORED, L_POWER as usize);
        let x = 19u32;

        let timer = Instant::now();
        let pairing_result = weil_pairing(&test_curve, &test_P, &test_Q, &new_order);
        println!(
            "first weil pairing elapsed : {} ms",
            timer.elapsed().as_millis()
        );
        let timer = Instant::now();
        let pairing_result2 = weil_pairing(
            &test_curve,
            &test_P,
            &test_curve.mul_small(&test_Q, x as u64),
            &new_order,
        );
        println!(
            "second weil pairing elapsed : {} ms",
            timer.elapsed().as_millis()
        );

        assert!(pairing_result.pow_big(&new_order).equals(&Fq::ONE) != 0);
        assert!(pairing_result2.pow_big(&new_order).equals(&Fq::ONE) != 0);

        assert!(pairing_result.pow_small(x).equals(&pairing_result2) != 0);
    }

    #[test]
    fn compute_l_power_pairing() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let basis_order = BigUint::from_slice(&BASIS_ORDER);
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let cofactor = &basis_order / &l_power;
        let (P, Q) = entangled_torsion_basis(&test_curve, &cofactor);
        assert!(test_curve.mul_big(&P, &l_power).isinfinity() != 0);
        assert!(test_curve.mul_big(&Q, &l_power).isinfinity() != 0);
        let ePQ = weil_pairing(&test_curve, &P, &Q, &l_power);
        println!("e(P, Q) = {}", ePQ);
        assert!(
            ePQ.pow_big(&l_power).equals(&Fq::ONE) != 0,
            "ePQ^d = {}",
            ePQ.pow_big(&l_power)
        );
        assert!(
            has_factored_order(ePQ, &[(2u32, L_POWER)]),
            "e(P, Q) has wrong order"
        );
    }
}
