use num_bigint::BigUint;

use crate::{
    fields::FpFESTAExt::BIGUINT_MODULUS,
    thetaFESTA::{Curve, Fq, Point},
};

/// Evalute the line g crossing the two points P and Q
/// and output g(R)
fn eval_line(E: &Curve, P: &Point, Q: &Point, R: &Point) -> Fq {
    let result;
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

    debug_assert!(V.isinfinity() != 0, "V is not identity");

    let q = BigUint::from_slice(&BIGUINT_MODULUS);
    let exp = (&q * &q - BigUint::from(1u32)) / order;
    f.pow_big(&exp)
}

#[cfg(test)]
mod tests {
    use crate::fields::FpFESTAExt::BASIS_ORDER;

    use super::*;
    use num_bigint::BigUint;
    use rand::prelude::*;
    use rand_chacha::ChaCha20Rng;
    #[test]
    fn compute_pairing() {
        let test_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let mut rng = ChaCha20Rng::from_entropy();
        let mut test_P = test_curve.rand_point(&mut rng);
        let mut test_Q = test_curve.rand_point(&mut rng);
        let basis_order = BigUint::from_slice(&BASIS_ORDER) * BigUint::from(2u32);
        let new_point = test_curve.mul_big(&test_P, &basis_order);

        assert!(new_point.isinfinity() != 0);

        // sample order for test
        let new_order = BigUint::from(41161u32);

        assert!((&basis_order % &new_order) == BigUint::from(0u32));

        let cofactor = &basis_order / &new_order;
        test_P = test_curve.mul_big(&test_P, &cofactor);
        test_Q = test_curve.mul_big(&test_Q, &cofactor);

        assert!(test_curve.mul_big(&test_P, &new_order).isinfinity() != 0);
        assert!(test_curve.mul_big(&test_Q, &new_order).isinfinity() != 0);

        let pairing_result = tate_pairing(&test_curve, &test_P, &test_Q, &new_order);

        assert!(pairing_result.pow_u64(41161u64, 16).equals(&Fq::ONE) != 0);

        println!("e(P, Q) : {}", pairing_result);
    }
}
