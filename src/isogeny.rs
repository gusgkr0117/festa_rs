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
                KummerLineIsogeny {
                    domain: *domain,
                    codomain: None,
                    kernel: *kernel,
                    degree,
                    edwards_multiples: vec![],
                }
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
                Aed = Aed.pow(&[self.degree as u8], 8) * prod_z;
                Ded = Ded.pow(&[self.degree as u8], 8) * prod_y;

                // Convert back into the Montgomery form
                A = Aed + Aed;
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

                self.codomain.unwrap().complete_pointX(&PointX::new_xz(&x_new, &z_new)).0
            }
        }

        /// Given a list of isogenies, evaluates the point for each isogeny in the list
        pub fn evaluate_isogeny_chain(P: &Point, isog_chain : &Vec<KummerLineIsogeny>) -> (Option<Curve>, Point) {
            if isog_chain.is_empty() {
                return (None, *P);
            }

            let mut Q = P.clone();
            for isog in isog_chain.iter() {
                Q = isog.evaluate_isogeny(&Q);
            }
            (Some(isog_chain.last().unwrap().get_codomain()), Q)
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
                        return vec![KummerLineIsogeny::new(E, Q, k)];
                    }

                    let mut k1 = (k * 8 + 5) / 10;
                    k1 = max(1, min(k - 1, k1));
                    
                    // Q1 <- l^k * Q
                    let mut Q1 = Q.clone();
                    for _ in 0..k1 {
                        Q1 = E.mul_small(&Q1, l as u64);
                    }

                    let mut L = recursive_sparse_isogeny(E, &Q1, l, k - k1);
                    let (coE, Q2) = evaluate_isogeny_chain(Q, &L);
                    let mut R = recursive_sparse_isogeny(&coE.unwrap(), &Q2, l, k1);

                    L.append(&mut R);
                    L
                }

                recursive_sparse_isogeny(start_curve, P, l, e)
            }

            let mut P = kernel.clone();
            let mut curve : Option<Curve>;
            let mut phi_list = Vec::new();

            for (l, e) in order.iter() {
                (curve, P) = evaluate_isogeny_chain(&P, &phi_list);
                curve = match curve {
                    Some(x) => Some(x), 
                    None => Some(domain.clone()),
                };
                let mut Q = P.clone();
                for (l1, e1) in order.iter() {
                    if l != l1 {
                        for _ in 0..*e1 {
                            Q = curve.unwrap().mul_small(&Q, *l1 as u64);
                        }
                    }
                }

                let mut psi_list = sparse_isogeny_prime_power(&curve.unwrap(), &Q, *l as usize, *e as usize);
                phi_list.append(&mut psi_list);
            }
            phi_list
        }
    };
}

pub(crate) use define_isogeny_structure;