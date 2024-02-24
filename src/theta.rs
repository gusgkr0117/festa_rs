#![allow(non_snake_case)]

// Macro for defining the following types:
// - ThetaPoint: An element of the level-2 theta structure, which is encoded by
// four projective coordinates: (Fq : Fq : Fq : Fq)
//
// - ThetaStructure: The parent of ThetaPoint, the identity point is the theta
// null point of type ThetaPoint. For arithmetic, this type also has an
// arithmetic precom. which can be reused for both doublings and isogenies.
//
// - product_isogeny: an implementation of an isogeny between elliptic products
// (E1 x E2) of type EllipticProduct = (Curve, Curve) with a kernel of type
// (CouplePoint, CouplePoint) where each CouplePoint represents a pair of points
// P1, P2 on E1 and E2 respectively

// Macro expectations:
// Fq      type of field element Fp^2
// Curve   type of elliptic curve in Montgomery model
// Point   type of point on a montgomery curve
// EllipticProduct    type of E1 x E2
// CouplePoint        type of point on E1 x E2
macro_rules! define_theta_structure {
    () => {
        use std::fmt;

        /// Given four elements of Fq, compute the hadamard transform using recursive
        /// addition.
        /// Cost: 8a
        #[inline(always)]
        fn to_hadamard(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
            let t1 = X + Y;
            let t2 = X - Y;
            let t3 = Z + T;
            let t4 = Z - T;

            let A = &t1 + &t3;
            let B = &t2 + &t4;
            let C = &t1 - &t3;
            let D = &t2 - &t4;
            (A, B, C, D)
        }

        /// Given four elements of Fq, first square each coordinate and
        /// then compute the hadamard transform
        /// Cost: 4S, 8a
        #[inline(always)]
        fn to_squared_theta(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
            let XX = X.square();
            let YY = Y.square();
            let ZZ = Z.square();
            let TT = T.square();

            to_hadamard(&XX, &YY, &ZZ, &TT)
        }

        // ========================================================
        // Functions for working with ThetaPoints
        // ========================================================

        /// Theta Point Struct
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPoint {
            X: Fq,
            Y: Fq,
            Z: Fq,
            T: Fq,
        }

        impl ThetaPoint {
            /// Use for initalisation, probably stupid, or at least should have
            /// a different name!
            pub const ZERO: Self = Self {
                X: Fq::ZERO,
                Y: Fq::ZERO,
                Z: Fq::ZERO,
                T: Fq::ZERO,
            };

            /// Compile time, create a new theta point from Fq elements
            pub const fn new(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                    T: *T,
                }
            }

            /// Create a new theta point from Fq elements
            pub fn from_coords(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                    T: *T,
                }
            }

            /// Recover the coordinates of the element
            pub fn coords(self) -> (Fq, Fq, Fq, Fq) {
                (self.X, self.Y, self.Z, self.T)
            }

            /// Recover the coordinates of the element
            pub fn list(self) -> [Fq; 4] {
                [self.X, self.Y, self.Z, self.T]
            }

            /// Compute the Hadamard transform of the point's coordinates
            pub fn hadamard(self) -> (Fq, Fq, Fq, Fq) {
                to_hadamard(&self.X, &self.Y, &self.Z, &self.T)
            }

            /// Square each of the point's coordinates and then
            /// compute the hadamard transform
            pub fn squared_theta(self) -> (Fq, Fq, Fq, Fq) {
                to_squared_theta(&self.X, &self.Y, &self.Z, &self.T)
            }
        }

        /// For debugging, pretty print out the coordinates of a point
        impl fmt::Display for ThetaPoint {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "{}\n{}\n{}\n{}\n", self.X, self.Y, self.Z, self.T)
            }
        }

        // ========================================================
        // Functions for working with ThetaStructures
        // ========================================================

        /// Theta Structure
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaStructure {
            null_point: ThetaPoint,
            arithmetic_precom: [Fq; 6],
        }

        impl ThetaStructure {
            /// Given the coordinates of a null point, create a null point and
            /// precompute 6 Fp2 elements which are used for doubling and isogeny
            /// computations.
            pub fn new_from_coords(X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> Self {
                let null_point = ThetaPoint::new(X, Z, U, V);
                Self {
                    null_point,
                    arithmetic_precom: ThetaStructure::precomputation(&null_point),
                }
            }

            /// Given a null point, store the null point and precompute 6 Fp2
            /// elements which are used for doubling and isogeny computations.
            pub fn new_from_point(null_point: &ThetaPoint) -> Self {
                Self {
                    null_point: *null_point,
                    arithmetic_precom: ThetaStructure::precomputation(null_point),
                }
            }

            /// Return the null point of the ThetaStructure
            pub fn null_point(self) -> ThetaPoint {
                self.null_point
            }

            /// For doubling and also computing isogenies, we need the following
            /// constants, which we can precompute once for each ThetaStructure.
            /// We require 6 Fq elements in total and the cost is
            /// 4S (sqr theta) + 1I + 15M (batch inversion) + 6M (calc)
            /// Cost: 1I + 21M + 4S
            #[inline]
            pub fn precomputation(O0: &ThetaPoint) -> [Fq; 6] {
                let (a, b, c, d) = O0.coords();
                let (AA, BB, CC, DD) = O0.squared_theta();

                // Use Montgomery's trick to match invert k values for a cost
                // of a single inversion and 3(k - 1) multiplications. Inversion
                // is done in place
                let mut inverses = [b, c, d, BB, CC, DD];
                Fq::batch_invert(&mut inverses);

                let y0 = &a * &inverses[0];
                let z0 = &a * &inverses[1];
                let t0 = &a * &inverses[2];

                let Y0 = &AA * &inverses[3];
                let Z0 = &AA * &inverses[4];
                let T0 = &AA * &inverses[5];

                [y0, z0, t0, Y0, Z0, T0]
            }

            /// Given a point P, compute it's double [2]P in place.
            /// Cost 8S + 6M
            #[inline(always)]
            pub fn set_double_self(self, P: &mut ThetaPoint) {
                let (mut xp, mut yp, mut zp, mut tp) = P.squared_theta();

                // Compute temp. coordinates, 8S and 3M
                xp = xp.square();
                yp = &self.arithmetic_precom[3] * &yp.square();
                zp = &self.arithmetic_precom[4] * &zp.square();
                tp = &self.arithmetic_precom[5] * &tp.square();

                // Compute the final coordinates, 3M
                let (X, mut Y, mut Z, mut T) = to_hadamard(&xp, &yp, &zp, &tp);
                Y *= &self.arithmetic_precom[0];
                Z *= &self.arithmetic_precom[1];
                T *= &self.arithmetic_precom[2];

                P.X = X;
                P.Y = Y;
                P.Z = Z;
                P.T = T;
            }

            /// Compute [2] * self
            #[inline]
            pub fn double_point(self, P: &ThetaPoint) -> ThetaPoint {
                let mut P2 = *P;
                self.set_double_self(&mut P2);
                P2
            }

            /// Compute [2^n] * self
            #[inline]
            pub fn double_iter(self, P: &ThetaPoint, n: usize) -> ThetaPoint {
                let mut R = *P;
                for _ in 0..n {
                    self.set_double_self(&mut R)
                }
                R
            }
        }

        // ========================================================
        // Compting the gluing (2,2)-isogeny from a product of
        // elliptic curves to a level 2 theta structure
        //
        // A lot of the code below is to compute a 4x4 matrix,
        // represented as an array [Fq; 16] to compute a symplectic
        // basis transformation to for the points into a form
        // compatible with the isogeny formula
        // ========================================================

        /// Given a point in the four torsion, compute the 2x2 matrix needed
        /// for the basis change
        /// M = [[a, b], [c, d]] represented as an array [a, b, c, d]
        /// Cost: 14M + 2S + 1I
        fn get_base_submatrix(E: &Curve, T: &Point) -> (Fq, Fq, Fq, Fq) {
            let (x, z) = T.to_xz();
            let (u, w) = E.x_dbl_coords(&x, &z); // Cost 3M 2S

            // Precompute some pieces
            let wx = &w * &x;
            let wz = &w * &z;
            let ux = &u * &x;
            let uz = &u * &z;
            let det = &wx - &uz;

            // Batch inversion
            let mut inverse = [det, z];
            Fq::batch_invert(&mut inverse);

            // Compute the matrix coefficients
            let d = &uz * &inverse[0]; // Computing d then a saves one negation
            let a = -&d;
            let b = -&(&wz * &inverse[0]);
            let c = &ux * &inverse[0] - &x * &inverse[1];

            (a, b, c, d)
        }

        /// Given the four torsion below the isogeny kernel, compute the
        /// compatible symplectic transform to allow the theta points to have
        /// a good representation for the gluing isogeny
        ///
        /// Input is expected to be K1 = (P1, P2), K2 = (Q1, Q2) in E1 x E2
        /// Inside (E1 x E2)[4].
        /// Cost 100M + 8S + 4I
        fn get_base_matrix(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
        ) -> [Fq; 16] {
            // First compute the submatrices from each point
            let (E1, E2) = E1E2.curves();
            let (P1, P2) = P1P2.points();
            let (Q1, Q2) = Q1Q2.points();

            // TODO: if these were submatrix computations were done together, we
            // could save 3 inversions... It would make the code harder to read
            // but would be an optimisation for the gluing.
            // Something to think about for when cost REALLY matters.
            // Cost: 4 x 14M + 2S + 1I = 56M + 8S + 4I
            let (g00_1, g01_1, g10_1, g11_1) = get_base_submatrix(&E1, &P1);
            let (g00_2, g01_2, g10_2, g11_2) = get_base_submatrix(&E2, &P2);
            let (h00_1, _, h10_1, _) = get_base_submatrix(&E1, &Q1);
            let (h00_2, h01_2, h10_2, h11_2) = get_base_submatrix(&E2, &Q2);

            // Compute the product of g1 * h1 and g2 * h2 as 2x2 matricies
            // and extract out the first column

            // first col of g1 * h1 = [[gh00_1, *], [gh10_1, *]]
            let gh00_1 = &g00_1 * &h00_1 + &g01_1 * &h10_1;
            let gh10_1 = &g10_1 * &h00_1 + &g11_1 * &h10_1;

            // first col of g2 * h2 = [[gh00_2, *], [gh10_2, *]]
            let gh00_2 = &g00_2 * &h00_2 + &g01_2 * &h10_2;
            let gh10_2 = &g10_2 * &h00_2 + &g11_2 * &h10_2;

            // start the trace with the identity
            let mut a = Fq::ONE;
            let mut b = Fq::ZERO;
            let mut c = Fq::ZERO;
            let mut d = Fq::ZERO;

            // T1
            a += &g00_1 * &g00_2;
            b += &g00_1 * &g10_2;
            c += &g10_1 * &g00_2;
            d += &g10_1 * &g10_2;

            // T2
            a += &h00_1 * &h00_2;
            b += &h00_1 * &h10_2;
            c += &h10_1 * &h00_2;
            d += &h10_1 * &h10_2;

            // T1+T2
            a += &gh00_1 * &gh00_2;
            b += &gh00_1 * &gh10_2;
            c += &gh10_1 * &gh00_2;
            d += &gh10_1 * &gh10_2;

            // Now we act by (0, Q2)
            let a1 = &h00_2 * &a + &h01_2 * &b;
            let b1 = &h10_2 * &a + &h11_2 * &b;
            let c1 = &h00_2 * &c + &h01_2 * &d;
            let d1 = &h10_2 * &c + &h11_2 * &d;

            // Now we act by (P1, 0)
            let a2 = &g00_1 * &a + &g01_1 * &c;
            let b2 = &g00_1 * &b + &g01_1 * &d;
            let c2 = &g10_1 * &a + &g11_1 * &c;
            let d2 = &g10_1 * &b + &g11_1 * &d;

            // Now we act by (P1, Q2)
            let a3 = &g00_1 * &a1 + &g01_1 * &c1;
            let b3 = &g00_1 * &b1 + &g01_1 * &d1;
            let c3 = &g10_1 * &a1 + &g11_1 * &c1;
            let d3 = &g10_1 * &b1 + &g11_1 * &d1;
            // 44M

            [a, b, c, d, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3]
        }

        /// Apply the base change described by M on a ThetaPoint in-place
        /// Cost: 16M
        #[inline]
        fn apply_base_change(P: &mut ThetaPoint, M: [Fq; 16]) {
            let (x, y, z, t) = P.coords();
            P.X = &M[0] * &x + &M[1] * &y + &M[2] * &z + &M[3] * &t;
            P.Y = &M[4] * &x + &M[5] * &y + &M[6] * &z + &M[7] * &t;
            P.Z = &M[8] * &x + &M[9] * &y + &M[10] * &z + &M[11] * &t;
            P.T = &M[12] * &x + &M[13] * &y + &M[14] * &z + &M[15] * &t;
        }

        /// Given a couple point as input, compute the corresponding ThetaPoint on
        /// the level two structure and then apply the basis change on this point
        /// Cost: 20M
        fn base_change_couple_point(P1P2: &CouplePoint, M: [Fq; 16]) -> ThetaPoint {
            let (P1, P2) = P1P2.points();
            let (mut X1, mut Z1) = P1.to_xz();
            let (mut X2, mut Z2) = P2.to_xz();

            // If we have the point (0, 0) swap to (1, 0)
            let P1_check = X1.iszero() & Z1.iszero();
            X1.set_cond(&Fq::ONE, P1_check);
            Z1.set_cond(&Fq::ZERO, P1_check);

            // If we have the point (0, 0) swap to (1, 0)
            let P2_check = X2.iszero() & Z2.iszero();
            X2.set_cond(&Fq::ONE, P2_check);
            Z2.set_cond(&Fq::ZERO, P2_check);

            // Take all products to get level-2 theta point
            let X = &X1 * &X2;
            let Y = &X1 * &Z2;
            let Z = &Z1 * &X2;
            let T = &Z1 * &Z2;
            let mut P = ThetaPoint::from_coords(&X, &Y, &Z, &T);

            // Finally apply the base change on the point
            apply_base_change(&mut P, M);
            P
        }

        /// For a theta point which is on an elliptic product,
        /// one of the dual coordinates will be zero. We need
        /// to identify the index of this zero element for the
        /// gluing isogeny codomain and evaluation functions
        fn zero_index(dual_coords: &[Fq; 4]) -> usize {
            let mut z_idx = 0;
            for (i, el) in dual_coords.iter().enumerate() {
                // When el is zero, the result is 0xFF...FF
                // and zero otherwise, so we can use this as
                // a mask for each step.
                let el_is_zero = el.iszero();
                z_idx |= (i as u32 & el_is_zero);
            }
            z_idx as usize
        }

        /// Given the 8-torsion above the kernel of order 2, computes the
        /// codomain ThetaStructure (2,2)-isogenous from a product of elliptic
        /// curves
        ///
        /// NOTE: this function is a little fussy as we need to avoid the
        /// zero-dual coordinate. There's a chance refactoring this could make
        /// it appear more friendly
        ///
        /// Cost: 8S 13M 1I
        fn gluing_codomain(
            T1_8: &ThetaPoint,
            T2_8: &ThetaPoint,
        ) -> (ThetaStructure, (Fq, Fq), usize) {
            // First construct the dual coordinates of the kernel and look
            // for the element which is zero
            // For convenience we pack this as an array instead of a tuple:
            let xAxByCyD: [Fq; 4] = T1_8.squared_theta().into();
            let zAtBzYtD: [Fq; 4] = T2_8.squared_theta().into();

            // One element for each array above will be zero. Identify this
            // element to get the right permutation below for filling arrays
            let z_idx = zero_index(&xAxByCyD);

            // Compute intermediate values for codomain
            let t1 = zAtBzYtD[1 ^ z_idx];
            let t2 = xAxByCyD[2 ^ z_idx];
            let t3 = zAtBzYtD[3 ^ z_idx];
            let t4 = xAxByCyD[3 ^ z_idx];

            // Invert all four values for codomain and images
            let mut inverse = [t1, t2, t3, t4];
            Fq::batch_invert(&mut inverse);

            // Codomain coefficients
            let mut ABCD = [Fq::ZERO; 4];
            ABCD[0 ^ z_idx] = Fq::ZERO;
            ABCD[1 ^ z_idx] = &t1 * &inverse[2];
            ABCD[2 ^ z_idx] = &t2 * &inverse[3];
            ABCD[3 ^ z_idx] = Fq::ONE;

            // Used for the image computation
            let a_inverse = &t3 * &inverse[0];
            let b_inverse = &t4 * &inverse[1];

            let (A, B, C, D) = to_hadamard(&ABCD[0], &ABCD[1], &ABCD[2], &ABCD[3]);
            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            (codomain, (a_inverse, b_inverse), z_idx)
        }

        /// Given a point and it's shifted value, compute the image of this
        /// point using the Fq elements precomputed during the codomain
        /// computation
        /// Note: The shift value is needed as one element of the dual coordinates
        /// is zero, so we take values from both and use linear algebra to recover
        /// the correct image.
        ///
        /// Cost: 8S + 10M + 1I
        fn gluing_image(
            T: &ThetaPoint,
            T_shift: &ThetaPoint,
            a_inv: &Fq,
            b_inv: &Fq,
            z_idx: usize,
        ) -> ThetaPoint {
            // Find dual coordinates of point to push through
            let AxByCzDt: [Fq; 4] = T.squared_theta().into();

            // We are in the case where at most one of A, B, C, D is
            // zero, so we need to account for this
            // To recover values, we use the translated point to get
            let AyBxCtDz: [Fq; 4] = T_shift.squared_theta().into();

            // We can always directly compute three elements
            let y = AxByCzDt[1 ^ z_idx] * a_inv;
            let z = AxByCzDt[2 ^ z_idx] * b_inv;
            let t = AxByCzDt[3 ^ z_idx];

            // To compute the `x` value, we need to compute a scalar, lambda,
            // which we use to normalise the given `xb`. We can always do this,
            // but we are required to compute an inverse which means checking
            // whether z is zero. If z is zero, we can compute lambda by simply
            // extracting the scaled zb and dividing by b from above. However,
            // when z is zero, we instead have to use t. To ensure that this is
            // constant time, we compute both using that inverting zero just
            // gives zero and conditionally swapping lanbda with lambda_t
            let zb = AyBxCtDz[3 ^ z_idx];
            let tb = &AyBxCtDz[2 ^ z_idx] * b_inv;

            let mut inverse = [zb, tb];
            Fq::batch_invert(&mut inverse);

            // Potentially one of these inverses are zero, but we do both
            // to avoid branching.
            let mut lam = &z * &inverse[0];
            let lam_t = &t * &inverse[1];
            lam.set_cond(&lam_t, z.iszero());

            // Finally we recover x
            let xb = AyBxCtDz[1 ^ z_idx] * a_inv;
            let x = xb * lam;

            // We now have values for `x,y,z,t` but to order them we need to use
            // the xor trick as above, so we pack them into an array with the
            // right ordering and then extract them back out
            let mut xyzt = [Fq::ZERO; 4];
            xyzt[0 ^ z_idx] = x;
            xyzt[1 ^ z_idx] = y;
            xyzt[2 ^ z_idx] = z;
            xyzt[3 ^ z_idx] = t;

            let (x, y, z, t) = to_hadamard(&xyzt[0], &xyzt[1], &xyzt[2], &xyzt[3]);

            ThetaPoint::from_coords(&x, &y, &z, &t)
        }

        /// Compute the gluing (2,2)-isogeny from a ThetaStructure computed
        /// from an elliptic product.
        fn gluing_isogeny(
            E1E2: &EllipticProduct,
            P1P2_8: &CouplePoint,
            Q1Q2_8: &CouplePoint,
            image_points: &[CouplePoint],
        ) -> (ThetaStructure, Vec<ThetaPoint>) {
            // First recover the four torsion below the 8 torsion
            let P1P2_4 = E1E2.double(&P1P2_8);
            let Q1Q2_4 = E1E2.double(&Q1Q2_8);

            // Use the four torsion to deterministically find basis change
            let M = get_base_matrix(&E1E2, &P1P2_4, &Q1Q2_4);

            // Take the points P1, P2 in E1 x E2 and represent them
            // as a theta point on a level 2 structure and map them
            // through the above basis change
            let T1_8 = base_change_couple_point(&P1P2_8, M);
            let T2_8 = base_change_couple_point(&Q1Q2_8, M);

            // Now it's time to compute the codomain and image of the isogeny
            // with kernel below T1, and T2.
            // We save the zero index, as we use it for the images, and we also
            // can precompute a few inverses to save time for evaluation.
            let (codomain, (a_inv, b_inv), z_idx) = gluing_codomain(&T1_8, &T2_8);

            // We now want to push through a set of points by evaluating each of them
            // under the action of this isogeny. As the domain is an elliptic product,
            // with elements of type CouplePoint, and the codomain is a ThetaStructure
            // with elements of type ThetaPoint, we need a new vector here and as we
            // iteratate through each CouplePoint, we can compute its image and push it
            // to the new vector.
            let mut theta_images: Vec<ThetaPoint> = Vec::new();

            // Per image cost =
            // 2 * (16M + 5S) for the CouplePoint addition
            // 2 * 20M for the base change
            // 8S + 4M + 1I for the gluing image
            // Total:
            // 76M + 18S + 1I per point
            for P in image_points.iter() {
                // Need affine coordinates here to do an add, if we didn't we
                // could use faster x-only... Something to think about but no
                // obvious solution.
                let P_sum_T = E1E2.add(P, &P1P2_4);

                // After we have added the points, we can use the gluing formula
                // to recover theta points on the level 2 theta structure. First we
                // must compute the basis change as we did for the kernel:
                let T = base_change_couple_point(&P, M);
                let T_shift = base_change_couple_point(&P_sum_T, M);

                // With a point and the shift value from the kernel, we can find
                // the image
                let T_image = gluing_image(&T, &T_shift, &a_inv, &b_inv, z_idx);
                theta_images.push(T_image);
            }

            (codomain, theta_images)
        }

        // ===================================================================
        // Compting general (2,2)-isogenies between theta structures
        //
        // NOTE: For the two steps before a product structure is reached, we
        // need additional symplectic transforms which is controlled by the
        // `hadamard` array of `bool`s. The purpose of these is to avoid null
        // points (or dual null points) which have zero elements, which are
        // incompatible with the doubling formula.
        // ===================================================================

        /// Given the 8-torsion above the kernel, compute the codomain of the
        /// (2,2)-isogeny and the image of all points in `image_points`
        /// Cost:
        /// Codomain: 8S + 13M + 1I
        /// Image: 4S + 3M
        fn two_isogeny(
            domain: &ThetaStructure,
            T1: &ThetaPoint,
            T2: &ThetaPoint,
            image_points: &mut [ThetaPoint],
            hadamard: [bool; 2],
        ) -> ThetaStructure {
            // Compute the squared theta transform of both elements
            // of the kernel
            let (xA, xB, _, _) = T1.squared_theta();
            let (zA, tB, zC, tD) = T2.squared_theta();

            // Batch invert denominators
            let mut inverse = [xA, zA, tB];
            Fq::batch_invert(&mut inverse);

            // Compute the codomain coordinates
            let mut A = Fq::ONE;
            let mut B = &xB * &inverse[0];
            let mut C = &zC * &inverse[1];
            let mut D = &tD * &inverse[2] * &B;

            // Inverses will be used for evaluation below
            let B_inv = &domain.arithmetic_precom[3] * &B;
            let C_inv = &domain.arithmetic_precom[4] * &C;
            let D_inv = &domain.arithmetic_precom[5] * &D;

            // Finish computing the codomain coordinates
            // For the penultimate case, we skip the hadamard transformation
            if hadamard[1] {
                (A, B, C, D) = to_hadamard(&A, &B, &C, &D);
            }
            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            // Now push through each point through the isogeny
            for P in image_points.iter_mut() {
                let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();
                if hadamard[0] {
                    (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                    (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
                } else {
                    (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
                }
                YY *= &B_inv;
                ZZ *= &C_inv;
                TT *= &D_inv;

                if hadamard[1] {
                    (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                }

                P.X = XX;
                P.Y = YY;
                P.Z = ZZ;
                P.T = TT;
            }

            codomain
        }

        /// Special function for the case when we are (2,2)-isogenous to a
        /// product of elliptic curves. Essentially the same as above, but with
        /// some small changes to deal with that a dual coordinate is now zero.
        /// Computes the codomain of the (2,2)-isogeny and the image of all
        /// points in `image_points`
        /// Cost:
        /// Codomain: 8S + 23M + 1I
        /// Image: 4S + 3M
        fn two_isogeny_to_product(
            T1: &ThetaPoint,
            T2: &ThetaPoint,
            image_points: &mut [ThetaPoint],
        ) -> ThetaStructure {
            // Compute the squared theta transform of both elements
            // of the kernel
            let (mut xA, mut xB, yC, yD) = T1.hadamard();
            (xA, xB, _, _) = to_squared_theta(&xA, &xB, &yC, &yD);

            let (mut zA, mut tB, mut zC, mut tD) = T2.hadamard();
            (zA, tB, zC, tD) = to_squared_theta(&zA, &tB, &zC, &tD);

            // Batch invert denominators
            let mut inverse = [xA, zA, tB, xB, zC, tD];
            Fq::batch_invert(&mut inverse);

            // Compute the codomain coordinates as well as precomputations for
            // the images
            let A = Fq::ONE;
            let B = &xB * &inverse[0];
            let C = &zC * &inverse[1];
            let D = &tD * &inverse[2] * &B;
            let B_inv = &xA * &inverse[3];
            let C_inv = &zA * &inverse[4];
            let D_inv = &tB * &inverse[5] * &B_inv;

            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            for P in image_points.iter_mut() {
                let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();

                (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);

                YY *= B_inv;
                ZZ *= C_inv;
                TT *= D_inv;

                P.X = XX;
                P.Y = YY;
                P.Z = ZZ;
                P.T = TT;
            }

            codomain
        }

        // ========================================================
        // Compting the symplectic transform to expose the
        // product structure and then compute the correct
        // splitting to Montgomery curves.
        // ========================================================

        // This function is a bit of a mess. Ultimately, we want to know whether
        // given some pair of indices whether we should multiply by minus one.
        // We do this by returning either 0: do nothing, or 0xFF...FF: negate
        // the value, which concretely is performed with set_negcond() on the
        // field element.
        //
        // Mathematically we have a few things to juggle. Firstly, although the
        // index should really be tuples (x, y) for x,y in {0,1} we simply index
        // from {0, ..., 3}. So there is first the identification of:
        //
        // 0 : (0, 0)
        // 1 : (1, 0)
        // 2 : (0, 1)
        // 3 : (1, 1)
        //
        // The next thing we need is the dot product of these indices
        // For example:
        // Take i . j is the dot product, so input (x, y) = (1, 3)
        // corresponds to computing:
        // (1, 0) . (1, 1) = 1*1 + 0*1 = 1
        //
        // This means evaluation of chi means the sign is dictated by
        // => (-1)^(i.j) = (-1)^1 = -1
        //
        // A similar thing is done for all pairs of indices below.
        //
        // TODO: there may be a nicer way to organise this function, but
        // I couldn't find a nice closed form for ±1 from a pair (i, j)
        // which i could compute on the fly without first matching from
        // x,y in {0,..,3} to i,j in {(0,0)...(1,1)} (which would mean
        // using a match anyway!!).
        fn chi_eval(x: &usize, y: &usize) -> u32 {
            match (x, y) {
                (0, 0) => 0,
                (0, 1) => 0,
                (0, 2) => 0,
                (0, 3) => 0,
                (1, 0) => 0,
                (1, 1) => u32::MAX,
                (1, 2) => 0,
                (1, 3) => u32::MAX,
                (2, 0) => 0,
                (2, 1) => 0,
                (2, 2) => u32::MAX,
                (2, 3) => u32::MAX,
                (3, 0) => 0,
                (3, 1) => u32::MAX,
                (3, 2) => u32::MAX,
                (3, 3) => 0,
                _ => 1,
            }
        }

        /// For a given index (chi, i) compute the level 2,2 constants (square).
        /// The purpose of this is to identify for which (chi, i) this constant
        /// is zero.
        fn level_22_constants_sqr(null_point: &ThetaPoint, chi: &usize, i: &usize) -> Fq {
            let mut U_constant = Fq::ZERO;
            let null_coords = null_point.list();

            for t in 0..4 {
                let mut U_it = &null_coords[t] * &null_coords[i ^ t];
                U_it.set_condneg(chi_eval(chi, &t));
                U_constant += &U_it;
            }
            U_constant
        }

        /// For each possible even index compute the level 2,2 constant. Return
        /// the even index for which this constant is zero. This only fails for
        /// bad input in which case the whole chain would fail. Evaluates all
        /// positions, and so should run in constant time.
        fn identify_even_index(null_point: &ThetaPoint) -> (usize, usize) {
            const EVEN_INDICIES: [(usize, usize); 10] = [
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
                (3, 0),
                (3, 3),
            ];
            // Initialise the return tuple
            let mut chi_zero = 0;
            let mut i_zero = 0;

            for (chi, i) in EVEN_INDICIES.iter() {
                let U_sqr = level_22_constants_sqr(null_point, chi, i);

                // When U_sqr is zero, U_sqr_is_zero = 0xFF...FF
                // and 0 otherwise, so we can use this as a mask
                // to select the non-zero index through the loop
                let U_sqr_is_zero = U_sqr.iszero();
                chi_zero |= (*chi as u32 & U_sqr_is_zero);
                i_zero |= (*i as u32 & U_sqr_is_zero);
            }
            (chi_zero as usize, i_zero as usize)
        }

        /// We can precompute 10 different symplectic transforms which
        /// correspond to each of the possible 10 even indicies which could be
        /// zero. We can select the right change of basis by using the above
        /// functions and then selecting the correct map accordingly.
        fn compute_splitting_matrix(null_point: &ThetaPoint) -> [Fq; 16] {
    #[rustfmt::skip]
            const MAPS: [[Fq; 16]; 10] = [
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::MINUS_ONE, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::MINUS_ONE, Fq::ZERO, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::ONE, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::ZETA,
                    Fq::ONE, Fq::MINUS_ZETA, Fq::ONE, Fq::ZETA,
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::MINUS_ZETA,
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::ZETA,
                ],
            ];

            // Identity the current location of the zero
            let zero_location = identify_even_index(null_point);

            // Compute the corresponding matrix to map the zero to
            // the desired place
            // TODO: is a match like this the best thing to do in Rust??
            let M: [Fq; 16];
            match zero_location {
                (0, 2) => M = MAPS[0],
                (3, 3) => M = MAPS[1],
                (0, 3) => M = MAPS[2],
                (2, 1) => M = MAPS[3],
                (0, 1) => M = MAPS[4],
                (1, 2) => M = MAPS[5],
                (2, 0) => M = MAPS[6],
                (3, 0) => M = MAPS[7],
                (1, 0) => M = MAPS[8],
                (0, 0) => M = MAPS[9],
                // The above locations are an exhaustive list of possible inputs, not sure how to tell rust this...
                _ => panic!("Unreachable"),
            }

            M
        }

        /// Map from a theta point to one which admits a splitting to elliptic
        /// products. Essentially requires computing the correct splitting
        /// matrix and then applying the isomorphism
        fn splitting_isomorphism(
            Th: ThetaStructure,
            image_points: &mut [ThetaPoint],
        ) -> ThetaStructure {
            // Compute the correct splitting matrix
            let mut O0 = Th.null_point();
            let M = compute_splitting_matrix(&O0);

            // Map the Theta Structure through the symplectic transform
            apply_base_change(&mut O0, M);

            // Map the points through the symplectic transform
            for P in image_points.iter_mut() {
                apply_base_change(P, M);
            }

            ThetaStructure::new_from_point(&mut O0)
        }

        /// Given a Theta point in the correct representation, compute two
        /// dimension 1 theta points.
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn split_theta_point(P: &ThetaPoint) -> ((Fq, Fq), (Fq, Fq)) {
            let (a, b, _, d) = P.coords();

            let P1 = (a, b);
            let P2 = (b, d);

            (P1, P2)
        }

        /// Given a dimension one null theta point, compute the corresponding
        /// elliptic curve in the Montgomery model by recovering the Montgomery
        /// coefficient A
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn null_point_to_montgomery_curve(O0: &(Fq, Fq)) -> Curve {
            let (a, b) = O0;

            let aa = a.square();
            let bb = b.square();

            let T1 = &aa + &bb;
            let T2 = &aa - &bb;

            let A = -&(&T1.square() + &T2.square()) / &(&T1 * &T2);

            Curve::new(&A)
        }

        /// Given a dimension one theta point, compute the corresponding
        /// elliptic curve point on the Kummer line (X : Z)
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn theta_point_to_montgomery_point(O0: &(Fq, Fq), P: &(Fq, Fq)) -> PointX {
            let (a, b) = O0;
            let (U, V) = P;

            let X = a * V + b * U;
            let Z = a * V - b * U;

            // TODO: rather than use this Type we could directly lift here instead...
            // exposing PointX to be public was the easiest way to keep everything from
            // eccore the same, which will help in the future
            PointX::new_xz(&X, &Z)
        }

        /// Given a ThetaStructure and set of points in the compatible form,
        /// compute the product of elliptic curves and affine points on these
        /// curves.
        fn split_to_product(
            Th: &ThetaStructure,
            image_points: &[ThetaPoint],
            num_image_points: usize,
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // First we take the domain theta null point and
            // split this to two level-1 theta null points
            let null_point = Th.null_point();
            let (O1, O2) = split_theta_point(&null_point);

            // Compute Montgomery curve from dimension one
            // null points
            let E3 = null_point_to_montgomery_curve(&O1);
            let E4 = null_point_to_montgomery_curve(&O2);
            let E3E4 = EllipticProduct::new(&E3, &E4);

            // Now compute points on E3 x E4
            let mut C: CouplePoint;
            let mut couple_points: Vec<CouplePoint> = vec![];

            for P in image_points.iter().take(num_image_points) {
                // Split to level 1
                let (P1, P2) = split_theta_point(P);
                // Compute the XPoint (X : Z) from each theta point
                let Q1X = theta_point_to_montgomery_point(&O1, &P1);
                let Q2X = theta_point_to_montgomery_point(&O2, &P2);

                // Lift these points to (X : Y : Z) on the curves
                let (Q1, _) = E3.complete_pointX(&Q1X);
                let (Q2, _) = E4.complete_pointX(&Q2X);

                // Package these points into a CouplePoint on
                // E3 x E4
                C = CouplePoint::new(&Q1, &Q2);

                // Push this into the output
                couple_points.push(C);
            }

            (E3E4, couple_points)
        }

        // ========================================================
        // Main Method! Compute the isogeny between elliptic
        // products
        // ========================================================

        // A general comment about "optimal strategies" -- For the isogeny we
        // have a kernel of order 2^n, and to make a step on the (2,2)-isogeny
        // graph what we want to do is scale this point to get 2^(n-1) * P,
        // which is a point of order two, which then allows us to compute the
        // codomain and images. However, as doubling is more expensive than
        // computing an image (in the general case) it's worth storing many
        // values along the way while doubling and pushing them all through the
        // isogeny. As these pushed through points have the same order, the
        // subsequent steps will need to be doubled less (at the cost of needing
        // to push through more points.)
        //
        // For a proper reference, see Sec 4.2 of
        // https://eprint.iacr.org/2011/506.pdf
        //
        // Gluing: Doubling cost: 16M 16S (We have to double two elliptic curve
        // points) Image cost: 76M + 18S + 1I (Very expensive!)
        //
        // All other steps: Doubling cost: 8S + 6M Image cost: 4S + 3M
        //
        // So ignoring the gluing step, we see images have 1/2 the cost
        // (mathematically this is expected as our doubling formula is
        // essentially just two isogeny images) and so the optimised strategy
        // is computed with a weight that doubling is 2X images.
        //
        // For a function to see how strategies are computed, see strategy.py
        // The current implementation "knows" that the gluing is more expensive
        // and so has extra costs for the leftmost branch of the tree.

        /// Compute an isogeny between elliptic products, naive method with no
        /// optimised strategy. Only here for benchmarking
        pub fn product_isogeny_no_strategy(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
            image_points: &[CouplePoint],
            n: usize,
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // Store the number of image points we wish to evaluate to
            // ensure we return them all from the points we push through
            let num_image_points = image_points.len();

            // Convert the &[...] to a vector so we can add points to this
            // dynamically during the optimal strategy
            let mut kernel_couple_pts = image_points.to_vec();

            // Include the kernel inside the vector of points
            // to evaluate. At each step, every element of the
            // vector should be evaluated
            kernel_couple_pts.push(*P1P2);
            kernel_couple_pts.push(*Q1Q2);

            // Compute points of order 8
            let P1P2_8 = E1E2.double_iter(&P1P2, n - 1);
            let Q1Q2_8 = E1E2.double_iter(&Q1Q2, n - 1);

            // Compute Gluing isogeny
            let (mut domain, mut kernel_pts) =
                gluing_isogeny(&E1E2, &P1P2_8, &Q1Q2_8, &kernel_couple_pts);

            // Do all remaining steps
            let mut Tp1: ThetaPoint;
            let mut Tp2: ThetaPoint;
            for k in 1..n {
                // Repeatedly double to obtain points in the 8-torsion below the kernel
                Tp1 = domain.double_iter(&kernel_pts[num_image_points], n - k - 1);
                Tp2 = domain.double_iter(&kernel_pts[num_image_points + 1], n - k - 1);

                // For the last two steps, we need to be careful because of the zero-null
                // coordinates appearing from the product structure. To avoid these, we
                // use the hadamard transform to avoid them,
                if k == (n - 2) {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, false])
                } else if k == (n - 1) {
                    domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
                } else {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, true])
                }
            }

            // Use a symplectic transform to first get the domain into a compatible form
            // for splitting
            domain = splitting_isomorphism(domain, &mut kernel_pts);

            // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
            // onto this product
            let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

            (product, couple_points)
        }

        /// Compute an isogeny between elliptic products, use an optimised
        /// strategy for all steps assuming doubling is always more expensive
        /// that images, which is not true for gluing.
        pub fn product_isogeny(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
            image_points: &[CouplePoint],
            n: usize,
            strategy: &[usize],
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // Store the number of image points we wish to evaluate to
            // ensure we return them all from the points we push through
            let num_image_points = image_points.len();

            // Convert the &[...] to a vector so we can add points to this
            // dynamically during the optimal strategy
            let mut kernel_couple_pts = image_points.to_vec();

            // Include the kernel inside the vector of points
            // to evaluate. At each step, every element of the
            // vector should be evaluated
            kernel_couple_pts.push(*P1P2);
            kernel_couple_pts.push(*Q1Q2);

            // Bookkeeping for optimised strategy
            let mut strat_idx = 0;
            let mut level: Vec<usize> = vec![0];
            let mut prev: usize = level.iter().sum();

            // =======================================================
            // Gluing Step
            // TODO:
            // Because of type differences there's annoying code reuse
            // for the optimal strategy here and again for every step
            // in the chain thereafter. Which is bothersome. Maybe there
            // is a better way to write this...
            // =======================================================
            let mut ker1 = *P1P2;
            let mut ker2 = *Q1Q2;

            while prev != (n - 1) {
                // Add the next strategy to the level
                level.push(strategy[strat_idx]);

                // Double the points according to the strategy
                ker1 = E1E2.double_iter(&ker1, strategy[strat_idx]);
                ker2 = E1E2.double_iter(&ker2, strategy[strat_idx]);

                // Add these points to the image points
                kernel_couple_pts.push(ker1);
                kernel_couple_pts.push(ker2);

                // Update the strategy bookkeepping
                prev += strategy[strat_idx];
                strat_idx += 1;
            }

            // Clear out the used kernel point and update level
            kernel_couple_pts.pop();
            kernel_couple_pts.pop();
            level.pop();

            // Compute Gluing isogeny
            let (mut domain, mut kernel_pts) =
                gluing_isogeny(&E1E2, &ker1, &ker2, &kernel_couple_pts);

            // ======================================================
            // All other steps
            // Compute the (2^n-1, 2^n-1)-chain in the theta model
            // =======================================================

            let mut Tp1: ThetaPoint;
            let mut Tp2: ThetaPoint;
            let mut kernel_len: usize;

            // Do all remaining steps
            for k in 1..n {
                prev = level.iter().sum();
                kernel_len = kernel_pts.len();

                Tp1 = kernel_pts[kernel_len - 2];
                Tp2 = kernel_pts[kernel_len - 1];

                while prev != (n - 1 - k) {
                    // Add the next strategy to the level
                    level.push(strategy[strat_idx]);

                    // Double the points according to the strategy
                    Tp1 = domain.double_iter(&Tp1, strategy[strat_idx]);
                    Tp2 = domain.double_iter(&Tp2, strategy[strat_idx]);

                    // Add these points to the image points
                    kernel_pts.push(Tp1);
                    kernel_pts.push(Tp2);

                    // Update the strategy bookkeepping
                    prev += strategy[strat_idx];
                    strat_idx += 1;
                }

                // Clear out the used kernel point and update level
                kernel_pts.pop();
                kernel_pts.pop();
                level.pop();

                // For the last two steps, we need to be careful because of the zero-null
                // coordinates appearing from the product structure. To avoid these, we
                // use the hadamard transform to avoid them,
                if k == (n - 2) {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, false])
                } else if k == (n - 1) {
                    domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
                } else {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, true])
                }
            }

            // Use a symplectic transform to first get the domain into a compatible form
            // for splitting
            domain = splitting_isomorphism(domain, &mut kernel_pts);

            // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
            // onto this product
            let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

            (product, couple_points)
        }
    };
} // End of macro: define_theta_structure

pub(crate) use define_theta_structure;

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::{
        ecFESTA::{CouplePoint, Curve, EllipticProduct, Fq, Point},
        fields::FpFESTAExt::{BASIS_ORDER, L_POWER, THETA_STRATEGY},
        pairing::weil_pairing,
        supersingular::point_has_factored_order,
        thetaFESTA::product_isogeny,
    };

    #[test]
    fn evaluate_theta_isogeny() {
        println!("Testing chain for FESTA parameters");

        let l_power = BigUint::from(2u32).pow(L_POWER);
        let new_order = &l_power * BigUint::from(4u32);
        let basis_order = BigUint::from_slice(&BASIS_ORDER);

        // Curve coefficients
        let A1_str = "d268d78b4a92566637452bf184ade847911397b523d84750d9beddf94d402fa5ca047b8a729844f8362094d212cd6fb9121a44241fbc8a0cc2bac1562744f015e3d41a6fa75014fa6a00ad4ac9ea34b202160d0876574324a65b7b2253142a6f408e89094e8ed1ba44fb56eae37cd6de728f177c83391c6763bbf493fc33874cee9259989d59d403886d48dad157f40238acb776ecf4e38fb5653ddada4afd6d0402099fc10c379b5190c8c195531815b4729da27d43cbaf41b5d53185cabaff515669eb38ef7a39eb1ecbae997f063f47fa9952aef62eed5ae8b647d578dcc1da25aaf197d1c2cfcb9043e9751c71fcf42ff71d5e975279fe6d0d1ee051b167d5e23b86f56b4cf99b5c52a8b22769b9b7c3d7d224afb3cb3e77d89837765a983c5cddc1e3174f55fbe955a3a4b13d87d9746edce0c31d1bac44eeca2cf0973c240c5b0c";
        let A2_str = "3d4cfe9cbfe6987fb8ea7a04f62c9ed28ef9f8ced54ed453efcdaa3621d00b77b0e7a7dc12b1d582032f0b98e174aa8dcbc514b88f39988a54e6ae6be1cf4f8746a86a8837fa6f10fb999732616484422cdd7ea7e501b858f6ac3f4cfdcfcd467301ec08817a9a2955e274391d16f57d0667eadbbbe77c7b48509ca21b5ccb254d879f13b10a2dd3d089aea44ecc8bb5efac9764f86630460a5a9ea718b4394dbb091ffa5bc28624f6004c4068a7b9dce26d84d02a5a0b50b1bc168d7b9ab63b920a1ab0f7b35c09e6ed717331281699648e4f723916033dc3784491c12470bec626c3a689b95b746c342cae065dd8a8592949148d0a4a38ce96238cb6c6a5bee6e272e67462a0afa95c4f06ad7151ae1742483d1ad3203718650d40907c717eed90dd7b4e3479b623b973cc1c981f57ea2c24e76a13c982a184a4252cc3e80d0d5a5e15";
        let (A1, _) = Fq::decode(&hex::decode(A1_str).unwrap());
        let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

        // Kernel Data
        let P1_X_str = "ab710145f609f21b96a9829e5e9ca411431f13fe7155e1a6424013967a326a691c701d3fe90cf1536cf504fa70e0c49fcc401b7e104bec8f5e648b64bf54694c83377ac40e795be89f394e2c581a79db88f9d58d5c3114d4c07488ed95467747ac88add4ddc2f8aa95eb47be54c115b294b34fa8fb8858224b5289cccb1b1fbef87cbb7b9f1acb012c9e4341b40e7db00bdb4ffcfe88a27175aa488e887c1cf1880cdaa0a1853d95790a0902169eee111a7cf76765d92f803f615a3fe73fc01d9d245bb0325dcb5f13017c344206b583d3324fd466509308f329c3e40a6ed152e01d763eeebc2525a3e1056f82831111b54684ebcf9776d4f0be11a92de06b5c47fca1188579177bca071a152e16133e17a055da98e6bc706557f35942c0e3f665a2db4c489ac91cbe484b3092fffba9b545799fcd59626fbfe5fb77c47aef6aeb5b7616";
        let P1_Y_str = "c05a766bc0b28b4977118f646c78d43a26c1f8b7ae75847f63f6ed94e4247d581e65e85a4d86eebbe8f7928670e8962e3c51cd843b7aec237c5e1d1a9a06e7a74c7a3672e0fd8f12704bb95ec930599e5dbb2790c9480f345feb74acb27050427107aba7227730e0a6f1951a3ff22e459ef26786d93778d6f464c47ff85773985726614a91f5da610a3cf5aad29f59bd70a7c3e41861eaa0e3636df16ec987f27816753310e1ad9014a3887b2e86055cc0b087fe71148fbef4ed961df36f7a74920904709290ac3215f35056257fb284c08b02516c6ba3f7c5b66001055f0f606fad2c24b4ab975012b10ef589e1600e90f02decefb45895aa28f70d931ccce1cb17a7a5383189e10927ec2ab1dc3fb159496c7f947fcd8f6159a04ae8fd1d401929978ebb0e7d2a7ea554cfb9230351d3b39f9da0975477e17ce9291f807fc39be24813";
        let P2_X_str = "495a0140eb8101aa90f13bdab5c65e53a655ef2b8d2280b14396f498dad675305b822928d54d6f81d494d5e639e28d6c29585e788454609dc69653e5c1294b6ffa69c17e5737a8bf506c0da45a792140d167e8d67fe0b39288f0cd67cb33989f147a8d88189967e384eb54d40b105424eedd368b89b8d26beac0ebbdb2fc2bb12aa105116ecd2d076fdea288f15af4d89f543a85077db767e2ca19c7c94f3331330d5d888cd541b6297558e74cb8da409ffd041094b0e8f6725fdcda3052dae881bf1f53bfab770946cc8167ae8dbd0c5352c4153a64f2e9999bd31ae1a8ce27ad9816d5dea6bb7101a856c3804e3d26eee3b5656c58ba43a00460a9fbd4f808dd21e998abe955513a11ee68e2bcd80769ed5a1876e2f7ddac97b50d790b494abd1aa733ff76956ed56fe2457c9175a8504a52bbc9ce1ac21b9a2c2b6fc517af15b59307";
        let P2_Y_str = "eeb2858418f1768adf8e61b803381587a72dd971a982a2cf40cc263905feddf6fd4131a134c89bcc46f8a758fb1ec03dfe4d973fbac5f2533774adbe8ef0be3b7efc4650923e36f79d906d443ec7fef2e97e1dd92189ebc0e2a29556dfa821508b353b48e8b231785ba7c95e38b9ad847b08e2949be98eacfdfaa3f7837b5c436dc3034d78b6d975bee83943981c25cabe5f33a6a39faea86c4501989cd51516f6033be07b5580b03020820669a0d84ab8456837f90c665526b41ec7fe76d3a4e0c8133890845b5b92532cfdd76c1171b2b824a9bbd7b5e1895798d313f2145abd944b91d7a0d4407b608584e72ac91e4d24472f670430bc7c7bf7c6e9829a4063a254d6dd245c52fe833494b430e4cc01f34a50d823a908427908f7ddd7704bbd0ad34bb7e3ae1d10061d214d6020ad2390730393ce8fd9da5afa8252923b88463ce20e";
        let Q1_X_str = "f1327865b3ecd10fd9af8a10a7301c031637d39d66da8caa9952cf5b9e1ab8b518303a2eb8dc913bef7fb352242e3d80b27056eedded2bc1fcf5e70bb5afd610ca066c1e34c64d8801d98af75c459f805a03719a7ac29ce07b5e5625a8eabab6c705b87d217fadee068dec01c9adf4071f5eace53cc7a50914af62dff0d5343c5efb28f67fa3f4018a7192eeb0ba11bcc034f7425bce1adb07369c92268a1f3d170547805c63ef9f44f099cf3f1733304fd012a3feec1ba33c7f9626e24c83ca089ff82005e39724dae7202ea80c5e3d1271b2aeceb6046dd26d01228385f8fb4493c05d0c0a6f9c3eb9fbe4ed0cd7dd2e9bf0896e6604b4a2f9e3b4ce8b7eb8330002f497bba8985854b8652b29a55b41af08ccc4e3db74d71d65ed7ecdef8ed821ae3ac65bfa60f23b089b674d21413a6bfc04cf0345ad6ee50d16a65b84e83564df0c";
        let Q1_Y_str = "6c3b1385160d8c41e4fcd4b0f12908a0ad58772a35b34fc0a15159a5a164ebca6f53873d8b44ea79435e1d3c7afa318286e2c952070f81cb46d8e55a9ce221fc0637e5f3e5dba2422e4858e3ce6720f5070e7ff93837c2d1df96327865497dfde7160e7e6fce952ddf6a4c355d0c65127840d73999e7b81a447d07e694d873639957fb67c8d3f5e10622a03e151899ec68f24c1879fa8f9e44f4749014e27990e60272fa04158db4ad5b7274eb88a456ed34d5e67cf1eaa8867f2158f85355c5ec68bb04b086622d28ef2ff62dfde9988d79e5444edadb59fbc999c9fdaa33488c94b178d9c218109682ee4ba51b8709fc81bb66f0683e219b6c5b9f2a7f05dda514329645d7858fa65093a86fe8864bf1dcabe9f5d3c952c27ab45a4491d9291f37b2876bb62e9d96ab3fa29bb5133c18ffab247e024cd30b23f2111eb7d4a3a9759b0a";
        let Q2_X_str = "34678391557c9e2776dff7177e4c6fd8e62442ecc0aceea3a2b7969b9720c493fb4d96affe30264e5ab4c5e73f63d35da8ccde60fe738470c54b849fa9a99d78a415e651c583f18990ce9483cc905c2338dd704b738c9158febea12db9bf2c546a101a0fbf6e4923f9fd89090ef4843a2005dab3c0a14f080b7589cec37f7ee9477da4cdadf3311ae9a599734ff2de0e9fb364627eb0f5ce10ec05cadd9f8e4fa904b646f943c1cd52a1ec46145364f1edc41700525db027b0014b97caab13f747fa72552e51e41e185064178d643ce49547eff8c06bb232ab372f6d4ecfeb8205ca7c42f70f0e8085d15a91d3801bdec8aee5ea1c516410df8fd3f0f09be42876614fc25ee675fc743d1d7f58b6a31af37370a5002ec639ed014f4a2c6fbf1b646b6a24b6cec49e3c9f7e1e1095e32f880640f47ed835c50aa7ea313d9068f1d7e6980a";
        let Q2_Y_str = "7c33493b0c0ab5bf88131c9bc9291b428d6812683d1f82b0f6c41dbc95e999b336a0a9efd0050255f0f290d76291ea7c72bb03d3e3130a2101ea5b4dea3a62397f56242e5ace05bf2754db2bd51ad90f37d92a4b7fdf9e705b4788b2a9e6088faf446475bfeb4cc8fd19b6655c5953611e673d3cd354132cf484e9e395eb557991ddd4d34a63036ac95801769d2ded726c8c77ae8d9bdc29dc42c34a75d03e0489047870e6f93ebcb5ec4cdf849b2c7f28404f3380c9b31b02e43e3c8a1c5e892aaa68d1c260f0d903055021027db1ac829ea2603f775e9c3b6d5594accb05c0656d9dcbd9c030510d1c3de017cce56e16d16a4e8a0aa5f31e1d9fc565cedd7757bac076f4d06fd9b24887af3d5829a0df783f0485890f9c1c1446cb99dedd1ddb81beb6f48f8c6a1928aa053d9c3c634ea104388d28a592e4c2aab38fd01999aa797c0b";
        let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
        let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
        let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
        let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
        let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
        let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
        let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
        let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

        // Points to push through isogeny
        let L11_X_str = "dea2bb5f803a0281b2c1dc948ce739aa77ccb66d167b62b069c4e6c5e8f2d6a21bb7b385d134f04608caf3f129aa29ed6ca8be9369ed27c43535bd282878c1e8c8907581110a22df7b2ab1358b8aa50fc665b7ca399aea781ef0ac250277583cb97dc7116f98e6a45810c468a9d5da4825f3af65c4a6b96debe9fa0b5f7805312eb5403a7fa97da650e475d1daa6a93db1e5bfef4ef4e977fae52f416a2af1cf2c04810a0c60deebb39bf19d01e74e0ed02cb85518d7285d804b56340e7e68bd3e59312c2b4cad5ece56415b5ec824c4be694df9255ff04cb3905af56adf4d6473efda1922ee615d1d3c0741eadf41d1c9587d49ceccaed909b8b0174e61ce6f0df628a03da0def573fd3b7d4850bfa751b73a37839ebb484f9a9f112d95d1ee2522dc3703c5102d957fc3dfbcc24d5618a6cf5b106ef1f23590fc603699aa5887f01302";
        let L11_Y_str = "abe10e6a88932185ab71f34c0a808b1a6e34844f6bea614c04515b2755f84947d13eb29658659a4d4c66a0e2031964448f69b14a3cfe8a08533981747a82d1e648bfd9faba0b5349fe1a057e356b4fc9fa6d71078173665cfb07e19b03e2da3c45aab7e70ca1c80a0a9e0bcb2b551103c58da65277897c2470b70a05595ade7f7950c04f57497c8421d9081591529f54942fd749fb86c28ebedfd85e5d3e8f92131275c689849076add344565b166beee4d889791d9072f313305ce5c8c1f39c6c339e986d852672a1eff6a7f08ed7c99abf884db79806adf520c2fc8a34d6ddffa600454f181b563375de64f62e16ac48aa413301d5d4443ad325552d3381c4b18e317ed659978967a71d3a68c9238b6828d4f9483f5c298a872203036e5380112b1da86f0f7b44c205d20c1290eec4c644f0b5aefa1dbaac61142bc8d4232ac10cc30c";
        let (L11_X, _) = Fq::decode(&hex::decode(L11_X_str).unwrap());
        let (L11_Y, _) = Fq::decode(&hex::decode(L11_Y_str).unwrap());

        let L12_X_str = "67bd25d443ab35439cde45a4cae9deede7c7ab808d3cf080ae2d2fa8ddfbbde33fb472b24d0838a80ee4af7f4d3fe5ae85a1becad033785fb7024b4d203b3372c9a71816f0f16108fe4f9f989dec5d8427af9f49a7f3c43c3f7e8eef91bd8126cdaf3e953ae90fcb679e9d84beca4b803af9f8f6ab0af82ab00693a258121450afa5fb59227855a0e6b7301362382b8bd040d72d5e48a5627bd73d1b8c4c3e94d00fe3807cf8708a8e5752c0c4b42e9b68cf6640bda42a05efb996560aac1c407a6c0b36bbbcb9c6355cb897bf7191aafdd8ab7582becd10a5f57e40851debd9076ad86d4c425152e8a5609ee2a2c50737d05e09a2089cec0066d6e8a8f8e613aa06cd268e09725ca0ea0ad865753779e2350f7324a41145f57eec8ca3527fcd28024cab1cf831ae25398f88ad3052cd6a8955573af4e79e77eea6f3353282234072140b";
        let L12_Y_str = "fc5fb5ec8552da8cdde3fb67034abcd9832e2ae0ee73a46de6d17d1f3565eeaffb053011a34e951572ac8ba83840f713c104f171fd76979ccdb3fdba531df2fe172d37bd01a7246a3132887c6ec5cc76a6e8a031df962a86bbbdfc0ec813b501a1246a34d784a4c4bca9cb10785a7eea9a068df98907769c351292287cfa6f441de196f2406e966ab412140789260f8c918390c8149df0513bf3fda0eb1d275dcc00d30f507ce8ad17a91cdb083227b971162b04284eee36506ccf6c71b1ee07d3f2ac06e74192cd560e3d6093a4375b18b7b6281c9e1cc0583495e5169e85858ac2be3e600c6a8ba53f6c7bd9ede931fca251fa0fe2a218ae3d47565d45fdcdbaea8a262a22f9b472ee2c6fb3df9963fa2a6d8899485659ea42110f0b5dc8c92855de8dee7e64ebae6c45af01f83e36a7a3f8775808867cb165e04918525889d1661706";
        let (L12_X, _) = Fq::decode(&hex::decode(L12_X_str).unwrap());
        let (L12_Y, _) = Fq::decode(&hex::decode(L12_Y_str).unwrap());

        let L21_X_str = "5b3bcf51421852af967115e3824bbeffe9b9aab16762e530e35c32a77b0b37c62a7e0ce96f7af6b831cf7193f807ed48749ebee12c5fcf8a2fd9bf187b5a58856888f99c81fb89644138d0b136589aa28da0feb085cf3d2a43c0d54bd8c87aec9b818542ff96c0a95643402026f29eba75c9cdf0deabd841bc1f1e88f5dd04bdac64174c5b7dc7af8428f3c74c828d5dd37ff1afda9e00d4a4b8dc3d579bb1425809d8a751619517a25f9ac408ddc3391c73b170127d1e2492d8877ad069052b6cba9159ff2590bc8e628496002bca3665174d414ed67d842f93931418560264ebb0359f9ae42ca64ddd75834cdeae6e7e7172f981a419406fe24401c4d70198cfadb6f301125b7c3f37a65babab3f347c935d0cbf8ef8e7dda642767960e73380fef72976e386007bba68fefc04ab876a50db3e0c7f9b12965b61e4fbde8f887c47f70d";
        let L21_Y_str = "5ff5ef87185ee6951c538ef7f5e16b17c3c9a28b8e7828743e45472482d01062d5d616cfb34d86fb0a66cf65aa7820353b7d8df025e91b6685009ef4b2231983c3b428e79e11afe00aa7f69a7b4f4f1c09f51953b44ba84b2b0d5e80c47f66a5e883353acca8b79cd65b5c77185c028a5af8b08cdc33b01ae35be1866b932cce8ef2803800a50f07dd6d835503102ed21871f13814051ae73d8d1e29ef85a94032082996a0c8e58c5460b9dfef65311511e8030bb434ce7e4bd98ee7981143f2ac175e141776e583ce6ad2e2e06bc6ad2e20effe7969b9cd57ac1adb7cd89decce7347fcf1522dfed1ff2663d3942fe9fd86dd798d9ff6edb0db17472259e5b63086e4e8e04f112c474f8a4b64dd5a89cfb176ea9742f954a507b21dc3163bc1609dd5077daec60c14851496d3e56612210a8bff825890752462c02e8b720d909b685308";
        let (L21_X, _) = Fq::decode(&hex::decode(L21_X_str).unwrap());
        let (L21_Y, _) = Fq::decode(&hex::decode(L21_Y_str).unwrap());

        let L22_X_str = "4eb3160cce9bb8a02fbd87d6f070cfdab2f0465f2ed692457f56a9699f3d2a51385d1b0ea6c112676da50db2cd916295aeae5a0c8fab327c9019ecdbbae67d0677c38688b7988a99bd07bc4f67f47de987f535bdf4f7af89b6f17aebd5841eebc675a3bd0622b4445a86f469ab6478760e8d54f4810e60f82f25df97957768708b2373e28c5bdc051f5eb829ff63e25f024f546d7c0d2d30a90d235961ff861f14163fe43d8b8c9b8fc5f8bc8ee7970f3aa23ab78eb4a9a836ff087806843a4935860aab038b0809f7df120f5eeaa5fcf6c936158448cac458ac9c7ca3e89d830efde790240a79a515c4cb76e6e93ed2b121bbb75397153154ed16345358a9b2e6e28c70cef892dc01fb13ee017d04ee3848cf67fc8107b80b02b139ac9cb987c47cc3b8746a49d9e21e68643861d94cfbdce8aaed2ec3fd81be6b4aeefd27a3d85d7c12";
        let L22_Y_str = "6eb74be4d694cfd39f3a216ea77514079c66b8dff598af6588683946493a691d812685e352c218f95f465e2b21da95012c409cbf5c9925ae5e2c3ba6249164641d0314f03de05f3f9cd6851ae2681f04a7d4d04a76999ff6ce0dcb73f2d1e76fe9c86a403c62f4de281c3bccbfa72410a9bf26a7d2959577be8d9b23de9c634fd096b468e1da9ee98a4a449ac99bd590c1f9b60615dc748b6d7966e61e373edb7b07c97cffcd09ebc8daa380977a10555f800229e2bd158acc311e4e0bda66dae5f2b6d0745d62ed3d9a7d5415ba6aecfcc7db6eeff54476179c59e5cdb6f62c32f634b2de573181629b25d3febea4943982f72b535389183118caf16685fdb5ccfa73e609f9a7a98b7b1209231d0d9004cdff6b425cf578a4f586473e9a4fbf32b6fbae0cd1ecc1ea48d6ec950b2c42fe9a841a9cd559ab0154dc60b73ea290d1d37613";
        let (L22_X, _) = Fq::decode(&hex::decode(L22_X_str).unwrap());
        let (L22_Y, _) = Fq::decode(&hex::decode(L22_Y_str).unwrap());

        // Result to check against (Debugging)
        let im_L11_X_str = "5d1e75b96024eb526a1a8ab9ac93c105b12107125965d5005557a5255077a56e4b435a66fc17d554cc6b0639811ede4e7a06f4f93bf2ff9858115305eb3895913a737109c7b3e446cc66a7fe4b6a0836cb46c8c39d542a6fcb7d7585c5753196b80f2e53c8f385ff95dbf841303818fd95545a9e7bf6216cf2fac12d544827310b65c9c729bfea6c8ae27b8ab5947e183656eb8391f0bff1b983c1f9f317227d17046a442f190320cf651d95b3fc3537dd5b21616e3218cd8706f861a35a8f09dab19205c574b4b3ecfd099238189ad3ab900928dd43641c4310a45715cc836d76d6863e54847af02a5dd9c29adfb067b0c7f16f8a2d2d35c64427b8606e62e59267938eb9e8d1047b3676aacae48d55af42f0fd5a250828b88333c7a4c87846b9a22d489d6deb3746ee7ae9d173cbb4bb66071fca218cbab6c64b3e2b6fe230efcb4e08";
        let im_L11_Y_str = "c25bd83bb502e6a8ffdbafadd577c2529b457cdf1d470b764b61f05a0a5a331c20642902f1709d48295d691e8f6b47b423e653ac2054b59c04f980aa07586b5e24213dd335e7262b91c4c1b620817802969238a1f459b7f9bdc2ce31e5a38d5fce542132dc234c1645d95f171e7aea4a216a11b9eea3ee655c4e35f71f74a77ffa84229b215f4e51b8381637cf04817225f6b2333aa5de78f8af4f8d757883847e04d77ece30e000a72050c4664fd01a7c8e0c8afe4a4c2cb3ee56fe187effc9401a275f88969f7649d30bfcd1f23b64521b7dce274cca02289c3e41935d87a0a3a52d36342d617f2ae63cfdc0cc6ae2c20f11a61da9b3c695609a2ac4f818e4d78766597a0a44fc8b8428eedd62839a434318f5edc46db5279f2fe5794b89523ea97b3b7b532ad5b5c93d948ac175ab2131996ebce4d591408daa2296fda1210229df01";
        let (im_L11_X, _) = Fq::decode(&hex::decode(im_L11_X_str).unwrap());
        let (im_L11_Y, _) = Fq::decode(&hex::decode(im_L11_Y_str).unwrap());

        let im_L12_X_str = "b2b1c1489ad8010659ed45acee631ecd4f1654dcd014822d91963f53829e03dab4016ca088d377d9782cb983e78888b3f6b2fafbf9895d114dfc115899b08a4accc4ed467e3a80e2016542764f03afc3c0425c729eea94cf15738d1c50b68bcfd753e6f19a5bbbaf63abbec2b4ce10d29a255f9374c8bf718cf122e5b336c35872045afa5ca7c5d00ba1fc7b288bb3e8ff7729c2a66c97612ac497068324c662870240182b54358a2f3144c770800d05a4184ea4760161dc75a6a7d7604536dff6829af7578e8eab2eaf439eff87c69fd9d10a0ed06f34aff46cf15de71e6beaecaf7449b2e68a9f94f49b2e4fe98cc12c37e445b02762ab4fd3e8f75ad4ba1026da08d16f989b2a9cfb6949961e45b389ecff7406fb259821c2cee000fb1bde62619ceab4d3c1f4943c03f24a22ea0123fab70118c8c2da9d7540848cf7137bb06aef15";
        let im_L12_Y_str = "aa9e26c4724b97d8ce6a65ba2dc53d635c0105c6d1a7f6be7f3838714b5cd52b8c8c6d5a492331fabc513e1f759dd7417dcf4fa44b3354963ba32c6a2e3c22d0bd03484d96a980a9390acbbfb3243577f628d37865dac091b13b8c8f42b264a4b0ad3c1a5b6ca892f13eb3750aacf234e8db035ff6cda25b524a62a47b84da3e734f72e8aff2370e60ed2b2749c6d77085ee83e843bcd018660ae5e7458f1fd98b021559ea8727e553a4fef2200a9cc30c70d9a5aee5d8c7c35e0f9107ff9cf964e8c2d38fe56f7db19d067fe70fe4bcc50ad32e7034c62094731bffff7207e73362cd7be2048cc31242052f9fc77898d48ff8e1fb289720727d1594e4897edbabf8dd10c3d1894cc779483c12b3633341b1d464c58a22d21dfbfd8f24ddd7614c95a8d99ec8b6fcb0442f51fc3930e8944d15da2a8f588c36f61bc8c3da4a2815d17e11";
        let (im_L12_X, _) = Fq::decode(&hex::decode(im_L12_X_str).unwrap());
        let (im_L12_Y, _) = Fq::decode(&hex::decode(im_L12_Y_str).unwrap());

        let im_L21_X_str = "7b0af2b62ad1418d0e6f45e4f5e1e23aba8b21e603c1d327ecf65292fd1361210388c532d4076e3134a94ca28b936726f10a34062018ecdd278af21dc5aee3fe5129ed7592228bfb64d17ba67bb3dc02deefe95ce32e143c6e899eec2bd5b860364c9ca83fb06290d7a1fced8070049da7a2f61a279d79a2a6f3e3fc86efca46a8291e6b263c4df0cc086fa3db09584f321db9ae9c63d343657d3c8e9e2cb3186f0985a583085a80ff7c987d49e547c1ce65e99384a62a5b600939dbecff800962f22a508bccd7c87ce8fc0184fc75106534cd5a55f352a57b3cd3e55141deb2fc141768f2bcb667178496a0d89fd77947cc0d3fa0541b109fad7f39387f3eeca9ae1cd8cbb70efe29f6913dedd30b2a136f3660be6ffe1de52ea9da31ac9d392b667c6ff9feb5bee4c1b6da31ea62b9d484c7e63593324fa98982e6f2b4940a79471210";
        let im_L21_Y_str = "34652801a908dd0b5567b615b787a45b865b6a8fa4a5aecfbc7b434e671204705ee2a84de0c8ec82b60cfa259b6cbdde8f76a777e0e347bdd7c39c0cde925e61834fb6982de4d0cb1fa5b1144989013ec274b4b59878a0475a6dcf7d82d3b3104e2e2e7a01d5d1d308b705601c38368d3ef2b447c9eaef6629283a34fe157fff1e0392b3553f4da2a39ab66401d197239d897ed81f4d327ff448436ed2383a72f7121fb6e2674ed32f0b245dd96f50b5090c7acbffd9132da9675fd09322d74db9ffb487d5a9d9a259519cbd5d4f8d5f6553258e2f66ba1e8aa7bf543c8ec1ee32a2038d861c26b908a2de5ade9c60de6cb91cdddba7eda3928b129dcebb1d1976e3fc1a505baaf0bee43bfe3ca6bd826bfe37166e1d777a7310b8c58e7e136d894e730b77a68db35d3e026ded76d7bd80f9ebccc89ff0361084a8a2fc74c6dd01ba2c15";
        let (im_L21_X, _) = Fq::decode(&hex::decode(im_L21_X_str).unwrap());
        let (im_L21_Y, _) = Fq::decode(&hex::decode(im_L21_Y_str).unwrap());

        let im_L22_X_str = "6f8c129e738d35e908803a47aac5f45e0b6bca7ad4b0c59325d5169a1151ef8273f38e90e648cd8de58487544b8f2b9ae73697ba60d7ba9220bb25feffebc4638f7a7a8ded1f25367568e96d8c1e2bb38d4a9a13bd0ab02785695696a7d63a11f89622d2b373cf8c2412129812b30a62d032ffbcc192ed9398e1616e61c05a186d6ee3b3280ffb839387e7f4e0c0180a468cef504a39709e7ca1ed0172d0dd3655120857fe033936490b7c3ab0181a8c6b718c9ba86d36002c18dada113c735e5f628a4c0644bea3f78b16c2a9c16a5cfcdd5859827b5ce2c2752915d6c73b58f709ea3e52b695adca6a9c36ca9e4dc655145c9cdd2357e2cf969a5ad8486869d8e9357a0465345c1992a0117414f1dfc7138c5538ee7f6da349a8cbca8e599dc24064b5c162220bea04582fbe528a1f0abe9bca9698305c1f63054239f4dbfeb490b002";
        let im_L22_Y_str = "dd0683896fcf1abb06040dd596d5c25cc993a1776cb2d99f56c7212f029b8f68feb90cb54b3340396c189d42e619c2af1f42cee9a3c6ca7d5a35e7958e45bd203c6d6e96438c0a887e2ac0ad50e17a0ee24bc028d908a7edfcb9e9b9ff4d03556f2163090740ea97eac5de8d568f64ec6e7fac381b148b9f7570cda0ac71b01b0b561db179ce3a5adc29b9f4d0ac40d110fa01d7586f852296d99fb91d6feeea110c2d1a8499b556a372fbf5a7dee4d49e264fe452ba0ce592196957b9aeb226e7d053c645bc2a5c86b4955fb9b73bdaf480556532a58ed756f28630124725f8ef1b7fbc80f835175ba48baa9d39ffe9fc6d534aaa3d1631834f13a78cc1b356850ee4be86609e8eccc6e0355aa7ea26ec1f65cb5e6bf7b84a6ca5d76e03e72e59744bb1137a2c7b4688f0eeffb28360430f32589aa3fafba83a43cbfccd31a9645ed610";
        let (im_L22_X, _) = Fq::decode(&hex::decode(im_L22_X_str).unwrap());
        let (im_L22_Y, _) = Fq::decode(&hex::decode(im_L22_Y_str).unwrap());

        // Curves which define elliptic product
        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);
        let E1E2 = EllipticProduct::new(&E1, &E2);

        // Kernel Points on E1 x E2
        let P1 = Point::new_xy(&P1_X, &P1_Y);
        let P2 = Point::new_xy(&P2_X, &P2_Y);
        let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
        let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
        let P1P2 = CouplePoint::new(&P1, &P2);
        let Q1Q2 = CouplePoint::new(&Q1, &Q2);

        assert!(E1.on_curve(&P1));
        assert!(E1.on_curve(&Q1));
        assert!(E2.on_curve(&P2));
        assert!(E2.on_curve(&Q2));

        assert!(E1.mul_big(&P1, &new_order).isinfinity() != 0);
        assert!(E1.mul_big(&Q1, &new_order).isinfinity() != 0);
        assert!(E2.mul_big(&P2, &new_order).isinfinity() != 0);
        assert!(E2.mul_big(&Q2, &new_order).isinfinity() != 0);

        assert!(point_has_factored_order(&E1, &P1, &[(2u32, L_POWER + 2)]));
        assert!(point_has_factored_order(&E1, &Q1, &[(2u32, L_POWER + 2)]));
        assert!(point_has_factored_order(&E2, &P2, &[(2u32, L_POWER + 2)]));
        assert!(point_has_factored_order(&E2, &Q2, &[(2u32, L_POWER + 2)]));

        let eP1Q1 = weil_pairing(&E1, &P1, &Q1, &new_order);
        let eP2Q2 = weil_pairing(&E2, &P2, &Q2, &new_order);
        assert_eq!(
            eP1Q1 * eP2Q2,
            Fq::ONE,
            "The kernel is not maximally isotropic"
        );

        // Points to push through isogeny
        let L11 = Point::new_xy(&L11_X, &L11_Y);
        let L12 = Point::new_xy(&L12_X, &L12_Y);
        let L21 = Point::new_xy(&L21_X, &L21_Y);
        let L22 = Point::new_xy(&L22_X, &L22_Y);
        let L1 = CouplePoint::new(&L11, &L12);
        let L2 = CouplePoint::new(&L21, &L22);

        let image_points = [L1, L2];

        // Points to compare against
        let im_L11 = Point::new_xy(&im_L11_X, &im_L11_Y);
        let im_L12 = Point::new_xy(&im_L12_X, &im_L12_Y);
        let im_L21 = Point::new_xy(&im_L21_X, &im_L21_Y);
        let im_L22 = Point::new_xy(&im_L22_X, &im_L22_Y);

        // Compute chain
        let (E3E4, images) = product_isogeny(
            &E1E2,
            &P1P2,
            &Q1Q2,
            &image_points,
            L_POWER as usize,
            &THETA_STRATEGY,
        );

        println!("E1: {}", E1);
        println!("E2: {}", E2);

        let (E3, E4) = E3E4.curves();
        println!("E3: {}", E3);
        println!("E4: {}", E4);

        let (P1, P2) = images[0].points();
        let (P3, P4) = images[1].points();

        let P1_check: bool = P1.equals(&im_L11) | P1.equals(&(-im_L11)) == u32::MAX;
        let P2_check: bool = P2.equals(&im_L12) | P2.equals(&(-im_L12)) == u32::MAX;
        let P3_check: bool = P3.equals(&im_L21) | P3.equals(&(-im_L21)) == u32::MAX;
        let P4_check: bool = P4.equals(&im_L22) | P4.equals(&(-im_L22)) == u32::MAX;

        println!("First image matches what is expected: {}", P1_check);
        println!("Second image matches what is expected: {}", P2_check);
        println!("Third image matches what is expected: {}", P3_check);
        println!("Fourth image matches what is expected: {}", P4_check);
    }
}
