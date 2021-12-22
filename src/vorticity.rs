use crate::incompressibility::{nabla_pi_W, RADIUS_H};
use nalgebra::Vector3;

const EPSILON: f32 = 0.001;

/// This calculates ω_i as it is written in the paper.
#[inline(always)]
fn omega(
    positions: &Vec<Vector3<f32>>,
    velocities: &Vec<Vector3<f32>>,
    neighbours: &Vec<usize>,
    i: usize,
) -> Vector3<f32> {
    let mut sum = Vector3::zeros();

    for j in neighbours.iter() {
        let v_ij = velocities[*j] - velocities[i];
        // let v_ij = velocities[i] - velocities[*j];

        // Need to negate, ∇_pi W = -∇_pj W
        sum += v_ij.cross(&-nabla_pi_W(positions, i, *j, RADIUS_H));
    }

    sum
}

/// This gives ∇|ω_i|
/// Takes in omega for optimization (so I don't have to calculate omega twice)
#[inline(always)]
fn nabla_omega(omega: Vector3<f32>) -> Vector3<f32> {
    let omega_norm = omega.norm();
    
    if omega_norm != 0.0 {
        let omega_dx = omega.x / omega_norm;
        let omega_dy = omega.y / omega_norm;
        let omega_dz = omega.z / omega_norm;
        return Vector3::new(omega_dx, omega_dy, omega_dz);
    }

    return Vector3::zeros();
}

/// This gives N, which is η / |η| where η = ∇|ω_i|
/// Takes in omega for optimization (so I don't have to calculate omega twice)
#[inline(always)]
fn corrective_force(omega: Vector3<f32>) -> Vector3<f32> {
    let nabla_omega = nabla_omega(omega);

    if nabla_omega.norm() != 0.0 {
        return nabla_omega.normalize();
    }

    Vector3::zeros()
}

/// Calculates the vorticity force for particle i.
/// Denoted f^vorticity_i in the paper.
#[allow(non_snake_case)]
#[inline(always)]
pub fn f_vorticity(
    positions: &Vec<Vector3<f32>>,
    velocities: &Vec<Vector3<f32>>,
    neighbours: &Vec<usize>,
    i: usize,
) -> Vector3<f32> {
    let omega = omega(positions, velocities, neighbours, i);
    let N = corrective_force(omega);

    EPSILON * (N.cross(&omega))
}
