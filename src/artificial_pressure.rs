use crate::incompressibility::{w_poly6, RADIUS_H};
use na::Vector3;

// This is the value chosen for "n" in s_corr in the paper
// They recommended a value of 4
pub const S_CORR_N: f32 = 4.0;

// This is the value chosen for "k" in s_corr in the paper
// They recommended a value of 0.1
//pub const S_CORR_K: f32 = 0.01;
pub const S_CORR_K: f32 = 0.01;

// This is the value chosen for "|Δq|" in s_corr in the paper
// They recommended a value of 0.1h ... 0.3h
pub const S_CORR_DQ: Vector3<f32> =
    Vector3::new(0.2 * RADIUS_H, 0.0, 0.0);
    //Vector3::new(0.001, 0.0, 0.0);

/// The recommended pressure correction value to improve surface tension.
#[inline(always)]
pub fn s_corr(r: Vector3<f32>, h: f32) -> f32 {
    -S_CORR_K * (w_poly6(r, h) / w_poly6(S_CORR_DQ, h)).powf(S_CORR_N)
}
