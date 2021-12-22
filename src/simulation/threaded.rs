use crate::incompressibility::*;
use crate::{f_vorticity, viscosity, Boundary, Particle};
use na::Vector3;
use octree::Octree;
use rayon::iter::*;

/// Multithreaded implementation of Stage 1.
#[inline(never)]
fn stage_1_threaded(particles: &mut Vec<Particle>, dt: f32) -> Vec<Vector3<f32>> {
    // Constant gravity vector
    const G: Vector3<f32> = Vector3::new(0.0, -9.81, 0.0);

    particles.par_iter_mut().for_each(|particle| {
        // Assume mass is 1
        // The only external force is gravity
        // v_i <= Δt * f_ext(x_i)
        particle.velocity += dt * G;
    });

    let x_stars = particles
        .par_iter()
        .map(|particle| {
            // x*_i <= x_i + Δt * v_i
            particle.position + dt * particle.velocity
        })
        .collect();

    x_stars
}

/// Multithreaded implementation of Stage 2.
/// This one uses an octree.
#[inline(never)]
pub fn stage_2_octree_threaded(x_stars: &Vec<Vector3<f32>>) -> Vec<Vec<usize>> {
    let points: Vec<[f64; 3]> = x_stars
        .clone()
        .into_par_iter()
        .map(|x_star| [x_star.x as f64, x_star.y as f64, x_star.z as f64])
        .collect();

    let mut octree = Octree::new(points);
    octree.build((x_stars.len() as f32).sqrt() as usize);

    let h = RADIUS_H as f64;

    let indices = x_stars
        .par_iter()
        .enumerate()
        .map(|(i, x_star)| {
            let point = [x_star.x as f64, x_star.y as f64, x_star.z as f64];
            let octree_iter = octree.search(point, h);
            let mut x_star_indices = Vec::with_capacity(octree_iter.size_hint().1.unwrap() - 1);
            for (index, _) in octree_iter {
                if index != i {
                    x_star_indices.push(index);
                }
            }
            x_star_indices
        })
        .collect();
    indices
}

/// Multithreaded implementation of Stage 2.
/// This version uses a naive O(n^2) approach
#[inline(never)]
pub fn stage_2_naive_threaded(x_stars: &Vec<Vector3<f32>>) -> Vec<Vec<usize>> {
    let indices: Vec<Vec<usize>> = x_stars
        .par_iter()
        .enumerate()
        .map(|(i, x_star_i)| {
            let x_star_i = *x_star_i;
            x_stars
                .iter()
                .enumerate()
                .filter_map(|(j, x_star_j)| {
                    if i != j && (x_star_i - *x_star_j).norm_squared() <= RADIUS_H_POW2 {
                        Some(j)
                    } else {
                        None
                    }
                })
                .collect()
        })
        .collect();

    indices
}

/// Multithreaded implementation of Stage 3.
#[inline(never)]
fn stage_3_threaded(
    predicted_positions: &mut Vec<Vector3<f32>>,
    neighbour_indices: &Vec<Vec<usize>>,
    boundary: Boundary,
) {
    // Calculate λ_i's
    let lambdas: Vec<f32> = predicted_positions
        .par_iter()
        .enumerate()
        .map(|(i, _)| {
            let neighbours = &neighbour_indices[i];
            let lambda = lambda(&predicted_positions, neighbours, i);
            lambda
        })
        .collect();

    // Calculate Δp_i's
    let delta_p_s: Vec<Vector3<f32>> = predicted_positions
        .par_iter()
        .enumerate()
        .map(|(i, _)| {
            let neighbours = &neighbour_indices[i];
            let delta_p = delta_p(&predicted_positions, neighbours, &lambdas, i);
            delta_p
        })
        .collect();

    predicted_positions
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, predicted_position)| {
            *predicted_position += delta_p_s[i];

            predicted_position.y = predicted_position.y.clamp(0.0, f32::MAX);
            predicted_position.x = predicted_position.x.clamp(-boundary.x, boundary.x);
            predicted_position.z = predicted_position.z.clamp(-boundary.z, boundary.z);
        });
}

/// Multithreaded implementation of Stage 3.
#[inline(never)]
fn stage_4_threaded(
    particles: &mut Vec<Particle>,
    predicted_positions: Vec<Vector3<f32>>,
    neighbour_indices: Vec<Vec<usize>>,
    dt: f32,
) {
    // v_i <= 1/Δt * (x*_i - x_i)
    let mut new_velocities: Vec<Vector3<f32>> = particles
        .par_iter()
        .enumerate()
        .map(|(i, particle)| (predicted_positions[i] - particle.position) / dt)
        .collect();

    if std::env::args().find(|x| x == "nvort").is_none() {
        // First we apply vorticity confinement
        // Because this is multithreaded we need to build a Vec of vorticities first
        let vorticity_vec: Vec<Vector3<f32>> = new_velocities
            .par_iter()
            .enumerate()
            .map(|(i, _)| {
                f_vorticity(
                    &predicted_positions,
                    &new_velocities,
                    &neighbour_indices[i],
                    i,
                )
            })
            .collect();

        // Now we can apply vorticity
        new_velocities
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, new_velocity)| {
                // F = ma, and m = 1 in this implementation
                // So F = a. Hence, we can directly update velocity like so:
                *new_velocity += dt * vorticity_vec[i];
            });
    }

    if std::env::args().find(|x| x == "nvisc").is_none() {
        // Now we can apply viscosity
        // Because this is multithreaded we need to build a Vec of viscosities first
        let viscosity_vec: Vec<Vector3<f32>> = new_velocities
            .par_iter()
            .enumerate()
            .map(|(i, _)| {
                viscosity(
                    &predicted_positions,
                    &new_velocities,
                    &neighbour_indices[i],
                    i,
                )
            })
            .collect();

        // Now we can apply vorticity
        new_velocities
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, new_velocity)| {
                *new_velocity = viscosity_vec[i];
            });
    }

    // Finally we can update position.
    // To do so, we'll just assign all newfound values into the original particles Vec
    particles
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, particle)| {
            particle.position = predicted_positions[i];
            particle.velocity = new_velocities[i];
        });
}

pub fn run_sim_threaded(
    particles: &mut Vec<Particle>,
    dt: f32,
    iterations: usize,
    boundary: Boundary,
    octree: bool,
) {
    let mut x_stars = stage_1_threaded(particles, dt);
    let neighbour_indices = if octree {
        stage_2_octree_threaded(&x_stars)
    } else {
        stage_2_naive_threaded(&x_stars)
    };

    for _ in 0..iterations {
        stage_3_threaded(&mut x_stars, &neighbour_indices, boundary);
    }

    stage_4_threaded(particles, x_stars, neighbour_indices, dt);

    if std::env::args().find(|x| x == "assert").is_some() {
        for particle in particles.iter() {
            assert!(!particle.position.x.is_nan());
            assert!(!particle.velocity.x.is_nan());
        }
    }
}
