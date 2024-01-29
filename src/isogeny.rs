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
            pub fn evaluate_isogeny(self, P: PointX) -> PointX {
                let (xp, zp) = P.get_xz();
                let (p_sum, p_diff) = (xp + zp, xp - zp);

                let (mut x_new, mut z_new) = (Fq::ONE, Fq::ONE);
                for (ey, ez) in self.edwards_multiples.iter() {
                    let diff_ez = p_diff * ez;
                    let sum_ey = ey * p_sum;
                    x_new *= diff_ez + sum_ey;
                    z_new *= diff_ez - sum_ey;
                }

                PointX::new_xz(&x_new, &z_new)
            }
        }

        pub struct KummerLineIsogenyChain {
            domain : Curve,
            codomain : Option<Curve>,
            kernel : Point,
            degree : Vec<(u32, u32)>,
        }

        impl KummerLineIsogenyChain {
            pub fn new(domain: &Curve, kernel: &Point, degree: &Vec<(u32, u32)>) -> Self {
                KummerLineIsogenyChain {
                    domain: *domain,
                    codomain: None,
                    kernel: *kernel,
                    degree: *degree,
                }
            }
        }
    };
}

pub(crate) use define_isogeny_structure;
