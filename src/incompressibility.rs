use crate::artificial_pressure::s_corr;
use na::Vector3;
use std::f32::consts::PI;

// This is simply "ρ_0" in the paper
pub const RHO_0: f32 = 1.0;

// This is simply "h" in the paper
pub const RADIUS_H: f32 = 2.0;
pub const RADIUS_H_POW2: f32 = RADIUS_H * RADIUS_H;

// This is simply "ε" in the paper
pub const EPSILON: f32 = 10.0;

/// The W-Poly6 SPH density function.
#[inline(always)]
pub fn w_poly6(r: Vector3<f32>, h: f32) -> f32 {
    // We should only be calling this for neighbouring particles
    let r_norm = r.norm();
    if 0.0 <= r_norm && r_norm <= h {
        let coeff = 315.0 / (64.0 * PI * h.powf(9.0));
        return coeff * (h.powf(2.0) - r.norm_squared()).powf(3.0);
    }
    0.0
}

/// The standard SPH density estimator, ρ_i, used in the paper implementation.
#[inline(always)]
fn rho_i(positions: &Vec<Vector3<f32>>, neighbours: &Vec<usize>, i: usize) -> f32 {
    // Sum over all neighbour particles
    let mut sum = 0.0;
    for j in neighbours.iter() {
        // Authors mentioned to use poly6 for density estimation
        sum += w_poly6(positions[i] - positions[*j], RADIUS_H);
    }

    sum
}

/// The density constraint, C_i, per particle.
#[inline(always)]
#[allow(non_snake_case)]
pub fn C(positions: &Vec<Vector3<f32>>, neighbours: &Vec<usize>, i: usize) -> f32 {
    (rho_i(positions, neighbours, i)) / RHO_0 - 1.0
}

/// This gives ∇_pk Ci.
#[allow(non_snake_case)]
fn nabla_pk_Ci(
    positions: &Vec<Vector3<f32>>,
    neighbours: &Vec<usize>,
    k: usize,
    i: usize,
) -> Vector3<f32> {
    // First branch of piecewise
    // Here k = i
    if k == i {
        let mut sum = Vector3::zeros();

        for j in neighbours.iter() {
            sum += nabla_pi_W(positions, k, *j, RADIUS_H);
        }

        return sum / RHO_0;
    }

    // Second branch of piecewise
    // Here k = j

    // ∇_pj W is just the negative of ∇_pi W
    // So we'd have: -(-∇_pi W)
    // Which is effectively: ∇_pi W
    return nabla_pi_W(positions, i, k, RADIUS_H) / RHO_0;
}

/// This gives λ_i.
pub fn lambda(positions: &Vec<Vector3<f32>>, neighbours: &Vec<usize>, i: usize) -> f32 {
    let num = C(positions, neighbours, i);
    let mut denom_sum = 0.0;

    for k in neighbours.iter() {
        denom_sum += nabla_pk_Ci(positions, neighbours, *k, i).norm_squared();
    }

    // Need to also add the case where k = i
    denom_sum += nabla_pk_Ci(positions, neighbours, i, i).norm_squared();

    -(num / (denom_sum + EPSILON))
}

/// This gives ∇W, which should be the gradient w.r.t. p_i
/// Also: ∇_pj W = -∇_pi W
#[allow(non_snake_case)]
pub fn nabla_pi_W(positions: &Vec<Vector3<f32>>, i: usize, j: usize, h: f32) -> Vector3<f32> {
    // To calculate norm only once
    let r = positions[i] - positions[j];
    let r_norm = r.norm();

    // We should only be calling this for neighbouring particles
    if 0.0 < r_norm && r_norm < h {
        // Using the W spiky function for gradients, as mentioned in the paper
        let coeff = 15.0 / (PI * h.powf(6.0));

        // Really is just the beginning of the chain rule and r's l2 norm
        let derivative_shared = coeff * 3.0 * (h - r_norm).powf(2.0) / r_norm;

        let dw_dp_i_x = -derivative_shared * (r.x);
        let dw_dp_i_y = -derivative_shared * (r.y);
        let dw_dp_i_z = -derivative_shared * (r.z);

        return Vector3::new(dw_dp_i_x, dw_dp_i_y, dw_dp_i_z);
    }

    Vector3::zeros()
}

/// This gets Δp_i, the current iteration's calculated partial position
pub fn delta_p(
    positions: &Vec<Vector3<f32>>,
    neighbours: &Vec<usize>,
    lambdas: &Vec<f32>,
    i: usize,
) -> Vector3<f32> {
    debug_assert_eq!(positions.len(), lambdas.len());

    let mut sum = Vector3::zeros();

    if std::env::args().find(|x| x == "npress").is_none() {
        for j in neighbours.iter() {
            let r = positions[i] - positions[*j];
            sum += (lambdas[i] + lambdas[*j] + s_corr(r, RADIUS_H))
                * nabla_pi_W(positions, i, *j, RADIUS_H);
        }
    } else {
        for j in neighbours.iter() {
            sum += (lambdas[i] + lambdas[*j]) * nabla_pi_W(positions, i, *j, RADIUS_H);
        }
    }

    sum / RHO_0
}
