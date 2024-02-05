//! Macro for defining following types:
//! - isogeny: an implementation of isogeny between elliptic curves of Montgomery form
//!
//! Macro expectations:
//! Fq      type of field element Fp^2
//! Curve   type of elliptic curve in Montgomery form
//! Point   type of point on a montgomery curve
#[allow(non_snake_case)]

macro_rules! define_isogeny_structure {
    () => {
        /// Computation of isogenies between Kummer lines using x-only formula by Costello-Hisil-Renes
        /// This is used to compute only short prime degree isogenies
        /// To compute the composite degree isogeny, use KummerLineIsogenyChain
        #[derive(Clone)]
        pub struct KummerLineIsogeny {
            domain: Curve,
            codomain: Option<Curve>,
            kernel: Point,
            degree: usize,
            edwards_multiples: Vec<(Fq, Fq)>,
        }

        impl KummerLineIsogeny {
            /// Given domain, kernel and degree, compute the corresponding codomain curve
            /// and Edwards multiples
            pub fn new(domain: &Curve, kernel: &Point, degree: usize) -> Self {
                debug_assert!(kernel.isinfinity() == 0, "The kernel point cannot be zero");
                debug_assert!(domain.mul_small(kernel, degree as u64).isinfinity() != 0, "The order of the kernel point doesn't match");
                let mut result = KummerLineIsogeny {
                    domain: *domain,
                    codomain: None,
                    kernel: *kernel,
                    degree,
                    edwards_multiples: vec![],
                };

                result.compute_codomain_constants();
                result
            }

            /// Output its codomain curve
            pub fn get_codomain(&self) -> Curve {
                self.codomain.unwrap()
            }

            /// Computes the multiples that are used in both codomain computation and
            /// isogeny evaluation. We compute them once and reuse them for every evaluation.
            fn precompute_edwards_multiples(&mut self, d: usize) {
                let mut iter_point = Point::INFINITY;
                for _ in 0..d {
                    self.domain.addto(&mut iter_point, &self.kernel);
                    let (KX, KZ): (Fq, Fq) = iter_point.to_xz();
                    self.edwards_multiples.push((KX - KZ, KX + KZ));
                }
            }

            /// When the isogeny degree is odd, we compute the codomain using the
            /// Meyer and Reith Twisted Edwards trick (https://ia.cr/2018/782)
            pub fn compute_codomain_constants(&mut self) {
                // Extract Montgomery constant
                let (mut A, mut C): (Fq, Fq) = (self.domain.get_constant(), Fq::ONE);

                // Compute and store the edwards multiples
                self.precompute_edwards_multiples((self.degree - 1) / 2);

                // Compute the twisted Edward curve parameters of the domain curve
                let mut Ded = C + C;
                let mut Aed = A + Ded;
                Ded = A - Ded;

                // Compute product of Edwards multiples
                let (mut prod_y, mut prod_z) = (Fq::ONE, Fq::ONE);
                for (ey, ez) in self.edwards_multiples.iter() {
                    prod_y *= ey;
                    prod_z *= ez;
                }

                // Compute prod_y^8 and prod_z^8
                (prod_y, prod_z) = (prod_y * prod_y, prod_z * prod_z);
                (prod_y, prod_z) = (prod_y * prod_y, prod_z * prod_z);
                (prod_y, prod_z) = (prod_y * prod_y, prod_z * prod_z);

                // Compute the twisted Edward curve parameters of the codomain curve
                Aed = Aed.pow_big(&BigUint::from(self.degree)) * prod_z;
                Ded = Ded.pow_big(&BigUint::from(self.degree)) * prod_y;

                // Convert back into the Montgomery form
                A = Aed + Ded;
                C = Aed - Ded;
                A = A + A;

                // Set the computed codomain curve
                self.codomain = Some(Curve::new(&(A / C)));
            }

            /// Evalute the isogeny using Costello-Hisil formula
            /// for evaluating an odd degree isogeny on the Kummer line point P
            pub fn evaluate_isogeny(&self, P: &Point) -> Point {
                let (xp, zp) = P.to_xz();
                let (p_sum, p_diff) = (xp + zp, xp - zp);

                let (mut x_new, mut z_new) = (Fq::ONE, Fq::ONE);
                for (ey, ez) in self.edwards_multiples.iter() {
                    let diff_ez = p_diff * ez;
                    let sum_ey = ey * p_sum;
                    x_new *= diff_ez + sum_ey;
                    z_new *= diff_ez - sum_ey;
                }

                x_new = x_new * x_new * xp;
                z_new = z_new * z_new * zp;

                self.codomain.expect("No codomain curve").complete_pointX(&PointX::new_xz(&x_new, &z_new)).0
            }
        }

        impl fmt::Display for KummerLineIsogeny {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                if f.alternate() {
                    write!(f, "domain {} -> codomain {}", self.domain.j_invariant(), self.codomain.expect("No codomain curve").j_invariant())
                }else {
                    write!(f, "domain {} -> codomain {}", self.domain.get_constant(), self.codomain.expect("No codomain curve").get_constant())
                }
            }
        }

        /// Given a list of isogenies, evaluates the point for each isogeny in the list
        /// (domain, P) -> isog_1 -> ... -> isog_n -> (codomain, Q)
        pub fn evaluate_isogeny_chain(domain: &Curve, P: &Point, isog_chain : &Vec<KummerLineIsogeny>) -> (Curve, Point) {
            if isog_chain.is_empty() {
                return (*domain, *P);
            }

            let mut Q = P.clone();
            for isog in isog_chain.iter() {
                Q = isog.evaluate_isogeny(&Q);
            }
            (isog_chain.last().unwrap().get_codomain(), Q)
        }


        /// Computes a composite degree isogeny
        /// Given a kernel point P of degree l1^e1 * ... * lt^et
        /// compute the chain of Kummer isogenies
        pub fn factored_kummer_isogeny(domain : &Curve, kernel : &Point, order : &[(u32, u32)]) -> Vec<KummerLineIsogeny> {
            /// Compute chain of isogenies quotienting out a point P of order l^e
            pub fn sparse_isogeny_prime_power(start_curve: &Curve, P : &Point, l : usize, e : usize) -> Vec<KummerLineIsogeny> {
                /// Compute a chain of isogenies recursively
                pub fn recursive_sparse_isogeny(E: &Curve, Q : &Point, l : usize, k : usize) -> Vec<KummerLineIsogeny> {
                    use std::cmp::{max, min};
                    if k == 1 {
                        debug_assert!(E.mul_big(&Q, &BigUint::from(l)).isinfinity() != 0, "Q order is not 1");
                        return vec![KummerLineIsogeny::new(E, Q, l)];
                    }

                    let mut k1 = (k * 8 + 5) / 10;
                    k1 = max(1, min(k - 1, k1));
                    
                    // Q1 <- l^k * Q
                    let mut Q1 = Q.clone();
                    for _ in 0..k1 {
                        Q1 = E.mul_small(&Q1, l as u64);
                    }

                    debug_assert!(E.mul_big(&Q1, &BigUint::from(l).pow((k - k1) as u32)).isinfinity() != 0, "Q1 order is not k-k1");

                    let mut L = recursive_sparse_isogeny(E, &Q1, l, k - k1);
                    let (coE, Q2) = evaluate_isogeny_chain(E, Q, &L);
                    debug_assert!(coE.mul_big(&Q2, &BigUint::from(l).pow(k as u32)).isinfinity() != 0, "Q2 order is not k1");
                    let mut R = recursive_sparse_isogeny(&coE, &Q2, l, k1);

                    L.append(&mut R);
                    L
                }

                debug_assert!(e != 0, "zero exponential is invalid");

                recursive_sparse_isogeny(start_curve, P, l, e)
            }

            let mut P = kernel.clone();
            let mut curve : Curve = domain.clone();
            let mut phi_list : Vec<KummerLineIsogeny> = Vec::new();

            for (l, e) in order.iter() {
                let mut Q = P.clone();
                for (l1, e1) in order.iter() {
                    if l != l1 {
                        for _ in 0..*e1 {
                            Q = curve.mul_small(&Q, *l1 as u64);
                        }
                    }
                }
                debug_assert!(Q.isinfinity() == 0, "The kernel point must be nonzero");
                debug_assert!(curve.mul_big(&Q, &BigUint::from(*l).pow(*e)).isinfinity() != 0, "The order of the point is not correct");

                let mut psi_list = sparse_isogeny_prime_power(&curve, &Q, *l as usize, *e as usize);
                (curve, P) = evaluate_isogeny_chain(&curve, &P, &psi_list);
                phi_list.append(&mut psi_list);
            }

            debug_assert!(P.isinfinity() != 0, "P is not evaluated correctly");

            phi_list
        }
    };
}

pub(crate) use define_isogeny_structure;

#[cfg(test)]
mod tests {
    use crate::{
        ecFESTA::{evaluate_isogeny_chain, factored_kummer_isogeny}, thetaFESTA::{Curve, Fq}
    };

    use crate::fields::FpFESTAExt::BASIS_ORDER;
    use num_bigint::BigUint;

    use rand::prelude::*;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn compute_isogeny() {
        let start_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let mut rng = ChaCha20Rng::from_entropy();
        let randP = start_curve.rand_point(&mut rng);
        let randQ = start_curve.rand_point(&mut rng);
        let basis_order = BigUint::from_slice(&BASIS_ORDER) * BigUint::from(2u32);

        let factored_order : [(u32, u32);3] = [(3023, 1), (3359, 1), (4409, 1)];
        let new_order = BigUint::from(3023u32) * BigUint::from(3359u32) * BigUint::from(4409u32);
        let cofactor = &basis_order / &new_order;

        let kernP = start_curve.mul_big(&randP, &cofactor);
        let kernQ = start_curve.mul_big(&randQ, &cofactor);

        assert!(start_curve.mul_big(&kernP, &new_order).isinfinity() != 0);

        let isog_chain = factored_kummer_isogeny(&start_curve, &kernP, &factored_order);
        let (middle_curve, evalQ) = evaluate_isogeny_chain(&start_curve, &kernQ, &isog_chain);

        for isog in isog_chain.iter() {
            println!("{:#}", isog);
        }

        let isog_chain2 = factored_kummer_isogeny(&middle_curve, &evalQ, &factored_order);

        for isog in isog_chain2.iter() {
            println!("{:#}", isog);
        }
    }
}