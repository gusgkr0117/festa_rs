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
        use crate::pairing::weil_pairing;
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
                debug_assert!(
                    domain.mul_small(kernel, degree as u64).isinfinity() != 0,
                    "The order of the kernel point doesn't match"
                );
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

                self.codomain
                    .expect("No codomain curve")
                    .complete_pointX(&PointX::new_xz(&x_new, &z_new))
                    .0
            }
        }

        impl fmt::Display for KummerLineIsogeny {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                if f.alternate() {
                    write!(
                        f,
                        "domain {} -> codomain {}",
                        self.domain.j_invariant(),
                        self.codomain.expect("No codomain curve").j_invariant()
                    )
                } else {
                    write!(
                        f,
                        "domain {} -> codomain {}",
                        self.domain.get_constant(),
                        self.codomain.expect("No codomain curve").get_constant()
                    )
                }
            }
        }

        /// Given a list of isogenies, evaluates the point for each isogeny in the list
        /// (domain, P) -> isog_1 -> ... -> isog_n -> (codomain, Q)
        pub fn evaluate_isogeny_chain(
            domain: &Curve,
            P: &Point,
            isog_chain: &Vec<KummerLineIsogeny>,
        ) -> (Curve, Point) {
            if isog_chain.is_empty() {
                return (*domain, *P);
            }

            let mut Q = P.clone();
            for isog in isog_chain.iter() {
                Q = isog.evaluate_isogeny(&Q);
            }
            (isog_chain.last().unwrap().get_codomain(), Q)
        }

        /// Given the d-torsion basis <P, Q> = E[n]
        /// compute the image of the torsion basis up to the overall sign : ±(phi(P), phi(Q))
        pub fn evaluate_isogeny_chain_for_basis(
            domain: &Curve,
            basis: &(Point, Point),
            isog_chain: &Vec<KummerLineIsogeny>,
            torsion_order: &BigUint,
            isogeny_degree: &BigUint,
        ) -> (Curve, (Point, Point)) {
            let (P, Q) = *basis;
            let ((codomain, imP), (_, mut imQ)) = (
                evaluate_isogeny_chain(domain, &P, isog_chain),
                evaluate_isogeny_chain(domain, &Q, isog_chain),
            );

            let pair_E0 = weil_pairing(domain, &P, &Q, torsion_order);
            let pair_E1 = weil_pairing(&codomain, &imP, &imQ, torsion_order);

            // if e(P, Q)^n != e(phi(P), phi(Q))
            if pair_E0.pow_big(isogeny_degree).equals(&pair_E1) == 0 {
                imQ = -imQ;
            }

            (codomain, (imP, imQ))
        }

        /// Computes a composite degree isogeny
        /// Given a kernel point P of degree l1^e1 * ... * lt^et
        /// compute the chain of Kummer isogenies
        pub fn factored_kummer_isogeny(
            domain: &Curve,
            kernel: &Point,
            order: &[(u32, u32)],
        ) -> Vec<KummerLineIsogeny> {
            /// Compute chain of isogenies quotienting out a point P of order l^e
            pub fn sparse_isogeny_prime_power(
                start_curve: &Curve,
                P: &Point,
                l: usize,
                e: usize,
            ) -> Vec<KummerLineIsogeny> {
                /// Compute a chain of isogenies recursively
                pub fn recursive_sparse_isogeny(
                    E: &Curve,
                    Q: &Point,
                    l: usize,
                    k: usize,
                ) -> Vec<KummerLineIsogeny> {
                    use std::cmp::{max, min};
                    if k == 1 {
                        //debug_assert!(E.mul_big(&Q, &BigUint::from(l)).isinfinity() != 0, "Q order is not 1");
                        return vec![KummerLineIsogeny::new(E, Q, l)];
                    }

                    let mut k1 = (k * 8 + 5) / 10;
                    k1 = max(1, min(k - 1, k1));

                    // Q1 <- l^k * Q
                    let mut Q1 = Q.clone();
                    for _ in 0..k1 {
                        Q1 = E.mul_small(&Q1, l as u64);
                    }

                    //debug_assert!(E.mul_big(&Q1, &BigUint::from(l).pow((k - k1) as u32)).isinfinity() != 0, "Q1 order is not k-k1");

                    let mut L = recursive_sparse_isogeny(E, &Q1, l, k - k1);
                    let (coE, Q2) = evaluate_isogeny_chain(E, Q, &L);
                    //debug_assert!(coE.mul_big(&Q2, &BigUint::from(l).pow(k as u32)).isinfinity() != 0, "Q2 order is not k1");
                    let mut R = recursive_sparse_isogeny(&coE, &Q2, l, k1);

                    L.append(&mut R);
                    L
                }

                //debug_assert!(e != 0, "zero exponential is invalid");

                recursive_sparse_isogeny(start_curve, P, l, e)
            }

            let mut P = kernel.clone();
            let mut curve: Curve = domain.clone();
            let mut phi_list: Vec<KummerLineIsogeny> = Vec::new();

            for (l, e) in order.iter() {
                let mut Q = P.clone();
                for (l1, e1) in order.iter() {
                    if l != l1 {
                        for _ in 0..*e1 {
                            Q = curve.mul_small(&Q, *l1 as u64);
                        }
                    }
                }
                //debug_assert!(Q.isinfinity() == 0, "The kernel point must be nonzero");
                //debug_assert!(curve.mul_big(&Q, &BigUint::from(*l).pow(*e)).isinfinity() != 0, "The order of the point is not correct, (l,e) = ({}, {}), Point: {}", l, e, Q);

                let mut psi_list = sparse_isogeny_prime_power(&curve, &Q, *l as usize, *e as usize);
                (curve, P) = evaluate_isogeny_chain(&curve, &P, &psi_list);
                phi_list.append(&mut psi_list);
            }

            //debug_assert!(P.isinfinity() != 0, "P is not evaluated correctly");

            phi_list
        }
    };
}

pub(crate) use define_isogeny_structure;

#[cfg(test)]
mod tests {
    use crate::{
        discrete_log::{bidlp, ph_dlp},
        ecFESTA::{
            evaluate_isogeny_chain, evaluate_isogeny_chain_for_basis, factored_kummer_isogeny,
            KummerLineIsogeny,
        },
        fields::FpFESTAExt::{D1, D1_FACTORED, L_POWER},
        isogeny,
        pairing::weil_pairing,
        supersingular::{
            entangled_torsion_basis, generate_point_order_factored_D, has_factored_order,
            point_has_factored_order, torsion_basis,
        },
        thetaFESTA::{Curve, Fq},
    };

    use crate::fields::FpFESTAExt::BASIS_ORDER;
    use num_bigint::BigUint;

    #[test]
    fn evaluate_isogeny_for_basis() {
        let start_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let basis_order = BigUint::from_slice(&BASIS_ORDER);
        let l_power = BigUint::from(2u32).pow(L_POWER);
        let isogeny_degree_factored = [(2309u32, 1u32), (631u32, 1u32)];
        let isogeny_degree = isogeny_degree_factored
            .iter()
            .fold(BigUint::from(1u32), |r, (l, e)| {
                r * BigUint::from(*l).pow(*e)
            });
        let cofactor = &basis_order / &l_power;
        let (P, Q) = entangled_torsion_basis(&start_curve, &cofactor);
        let kernel = generate_point_order_factored_D(
            &start_curve,
            &isogeny_degree_factored,
            L_POWER as usize,
        );
        let phi = factored_kummer_isogeny(&start_curve, &kernel, &isogeny_degree_factored);
        let ((codomain, imP), (_, mut imQ)) = (
            evaluate_isogeny_chain(&start_curve, &P, &phi),
            evaluate_isogeny_chain(&start_curve, &Q, &phi),
        );

        let pair_E0 = weil_pairing(&start_curve, &P, &Q, &l_power);
        let pair_E1 = weil_pairing(&codomain, &imP, &imQ, &l_power);
        let pair_E0_d = pair_E0.pow_big(&isogeny_degree);

        println!("pair_E0^lb = {}", pair_E0.pow_big(&l_power));
        assert!(pair_E0.pow_big(&l_power).equals(&Fq::ONE) != 0);
        assert!(pair_E1.pow_big(&l_power).equals(&Fq::ONE) != 0);
        assert!(has_factored_order(pair_E0, &[(2, L_POWER)]));
        assert!(has_factored_order(pair_E1, &[(2, L_POWER)]));

        println!("pair_E0^d = {}", pair_E0_d);
        println!("pair_E1 = {}", pair_E1);
    }

    #[test]
    fn compute_isogeny() {
        let start_curve = Curve::new(&(Fq::TWO + Fq::FOUR));
        let basis_order = BigUint::from_slice(&BASIS_ORDER);

        let factored_D = [(2729u32, 1u32)];
        let factored_order1 = [(2309u32, 1u32)];
        let factored_order2 = [(631u32, 1u32)];

        let factored_whole = [(2309u32, 1u32), (631u32, 1u32)];
        let reversed_factored_whole = [(631u32, 1u32), (2309u32, 1u32)];

        let new_order1 = factored_order1
            .iter()
            .fold(BigUint::from(1u32), |r, (l, e)| {
                r * BigUint::from(*l).pow(*e)
            });
        let new_order2 = factored_order2
            .iter()
            .fold(BigUint::from(1u32), |r, (l, e)| {
                r * BigUint::from(*l).pow(*e)
            });
        let total_order = factored_whole
            .iter()
            .fold(BigUint::from(1u32), |r, (l, e)| {
                r * BigUint::from(*l).pow(*e)
            });

        assert!(
            (&basis_order % &new_order1) == BigUint::from(0u32),
            "The new order is wrong"
        );
        assert!(
            (&basis_order % &new_order2) == BigUint::from(0u32),
            "The new order is wrong"
        );

        let (P, Q) = torsion_basis(&start_curve, &factored_whole, L_POWER as usize);
        let R = generate_point_order_factored_D(&start_curve, &factored_D, L_POWER as usize);

        let middle_kernel = start_curve.mul_big(&P, &new_order2);
        let isog_chain1 = factored_kummer_isogeny(&start_curve, &middle_kernel, &factored_order1);
        let (middle_curve, imQ1) = evaluate_isogeny_chain(&start_curve, &Q, &isog_chain1);
        let (_, imP) = evaluate_isogeny_chain(&start_curve, &P, &isog_chain1);

        assert!(middle_curve.mul_big(&imP, &new_order2).isinfinity() != 0);
        println!("@0 : {}", middle_curve);

        let isog_chain2 = factored_kummer_isogeny(&middle_curve, &imP, &factored_order2);
        let (final_curve, imQ2) = evaluate_isogeny_chain(&middle_curve, &imQ1, &isog_chain2);
        let (_, imQ2) = evaluate_isogeny_chain(&start_curve, &imQ1, &isog_chain2);

        assert!(
            final_curve.on_curve(&imQ2),
            "phi(R) is not on the codomain curve"
        );
        println!("@1 : {}", final_curve);

        let isog_whole = factored_kummer_isogeny(&start_curve, &P, &factored_whole);
        let (_, im_whole_Q) = evaluate_isogeny_chain(&start_curve, &Q, &isog_whole);

        for (i, isog) in isog_whole.iter().enumerate() {
            println!("#{i} : {}", isog.get_codomain());
        }

        let isog_whole_dual =
            factored_kummer_isogeny(&final_curve, &im_whole_Q, &reversed_factored_whole);

        for (i, isog) in isog_whole_dual.iter().enumerate() {
            println!("#{} : {}", i + 2, isog.get_codomain());
        }
    }
}
