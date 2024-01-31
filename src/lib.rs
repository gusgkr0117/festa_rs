#![allow(non_snake_case)]

pub mod eccore;
pub mod festa;
pub mod fields;
pub mod finitefield;
pub mod isogeny;
pub mod pairing;
pub mod theta;

pub mod ecFESTA {
    pub type Fq = crate::fields::FpFESTAExt::Fp2;
    crate::eccore::define_ec_core! {}
    crate::isogeny::define_isogeny_structure! {}
}

pub mod thetaFESTA {
    pub use crate::ecFESTA::{
        CouplePoint, Curve, EllipticProduct, Fq, KummerLineIsogeny, Point, PointX,
    };
    crate::theta::define_theta_structure! {}
}
