use num_bigint::BigUint;

use crate::thetaFESTA::{Curve, Fq, Point};

/// Evalute the line g crossing the two points P and Q
/// and output g(R)
fn eval_line(E: &Curve, P: &Point, Q: &Point, R: &Point) -> Fq {
    let result;
    // When P or Q is the point at infinity i.e. D(g) = 0,
    // the dividing by this evaluation of the line must be ignored
    if P.isinfinity() != 0 || Q.isinfinity() != 0 {
        return Fq::ONE;
    }
    let (Px, Py) = P.to_xy();
    let (Qx, Qy) = Q.to_xy();
    let (Rx, Ry) = R.to_xy();
    let a = E.get_constant();
    if P.equals(&Q) != 0 {
        let lambda = (Fq::THREE * Px * Px + Fq::TWO * a * Px + Fq::ONE) / (Fq::TWO * Py);
        result = -lambda * Rx + Ry + (-Py + lambda * Px);
    } else if P.equals(&(-Q)) != 0 {
        result = Rx + (-Px);
    } else {
        let lambda = (Qy - Py) / (Qx - Px);
        result = -lambda * Rx + Ry + (-Py + lambda * Px);
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
        f = f * f * eval_line(E, &V, &V, Q) / eval_line(E, &E.double(&V), &E.double(&(-V)), Q);
        E.double_self(&mut V);
        if order.bit(t) {
            f = f * eval_line(E, &V, P, Q) / eval_line(E, &E.add(&V, P), &-E.add(&V, P), Q);
            E.addto(&mut V, P);
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
    use crate::{
        fields::FpFESTAExt::{BASIS_ORDER, L_POWER},
        supersingular::torsion_basis,
    };

    use super::*;
    use num_bigint::BigUint;
    #[test]
    fn compute_pairing() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let basis_order = BigUint::from_slice(&BASIS_ORDER);

        // sample order for test
        let new_order = BigUint::from(2u32).pow(L_POWER);

        assert!((&basis_order % &new_order) == BigUint::from(0u32));

        let (test_P, test_Q) = torsion_basis(&test_curve, &[(2u32, L_POWER)], 1);
        let x = 19u32;

        let pairing_result = weil_pairing(&test_curve, &test_P, &test_Q, &new_order);
        let pairing_result2 = weil_pairing(
            &test_curve,
            &test_P,
            &test_curve.mul_small(&test_Q, x as u64),
            &new_order,
        );

        println!("e(P, Q)^x : {}", pairing_result.pow_small(x));
        println!("e(P, nQ) : {}", pairing_result2);
        println!("e(P, Q)^r : {}", pairing_result.pow_big(&new_order));

        assert!(pairing_result.pow_big(&new_order).equals(&Fq::ONE) != 0);
        assert!(pairing_result2.pow_big(&new_order).equals(&Fq::ONE) != 0);

        assert!(pairing_result.pow_small(x).equals(&pairing_result2) != 0);
    }
}
