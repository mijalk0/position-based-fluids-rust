use crate::incompressibility::{w_poly6, RADIUS_H};
use na::Vector3;

// XSPH viscosity factor, denoted c in the paper.
// Authors recommended a value of 0.01
const C: f32 = 0.01;

/// Generates the viscosity update to velocity
#[inline(always)]
pub fn viscosity(
    positions: &Vec<Vector3<f32>>,
    velocities: &Vec<Vector3<f32>>,
    neighbours: &Vec<usize>,
    i: usize,
) -> Vector3<f32> {
    let v_i = velocities[i];

    let mut sum = Vector3::zeros();

    for j in neighbours.iter() {
        let v_ij = velocities[*j] - velocities[i];
        //let v_ij = velocities[i] - velocities[*j];
        let r = positions[i] - positions[*j];

        sum += v_ij * w_poly6(r, RADIUS_H);
    }

    v_i + C * sum
}
