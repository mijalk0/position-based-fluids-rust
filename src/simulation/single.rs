use crate::incompressibility::*;
use crate::{f_vorticity, viscosity, Boundary, Particle};
use na::Vector3;
use octree::Octree;

/// Stage 1:
/// This stage performs two updates: per particle.
/// 1) It applies the external force for a particle: v_i <= Δt * f_ext(x_i)
/// 2) It calculates (and returns) predicted positions: x*_i <= x_i + Δt * v_i
#[inline(never)]
fn stage_1(particles: &mut Vec<Particle>, dt: f32) -> Vec<Vector3<f32>> {
    let mut x_stars = Vec::with_capacity(particles.len());

    // Constant gravity vector
    const G: Vector3<f32> = Vector3::new(0.0, -9.81, 0.0);

    for i in 0..particles.len() {
        // Assume mass is 1
        // The only external force is gravity
        // v_i <= Δt * f_ext(x_i)
        particles[i].velocity += dt * G;

        // x*_i <= x_i + Δt * v_i
        let x_star_i = particles[i].position + dt * particles[i].velocity;
        x_stars.push(x_star_i);
    }

    x_stars
}

/// Stage 2:
/// This stage simply finds all neighbouring particles, using N_i(x*_i):
/// It returns a Vec, where Vec[i] is a vector of indices for neighbouring particles of particle[i]
/// This version uses an octree
#[inline(never)]
pub fn stage_2_octree(x_stars: &Vec<Vector3<f32>>) -> Vec<Vec<usize>> {
    let mut indices = Vec::with_capacity(x_stars.len());

    let points: Vec<[f64; 3]> = x_stars
        .clone()
        .into_iter()
        .map(|x_star| [x_star.x as f64, x_star.y as f64, x_star.z as f64])
        .collect();

    let mut octree = Octree::new(points);
    octree.build((x_stars.len() as f32).sqrt() as usize);

    let h = RADIUS_H as f64;

    for i in 0..x_stars.len() {
        let mut x_star_indices = Vec::new();
        let point = [
            x_stars[i].x as f64,
            x_stars[i].y as f64,
            x_stars[i].z as f64,
        ];
        for (index, _) in octree.search(point, h) {
            if index != i {
                x_star_indices.push(index);
            }
        }
        indices.push(x_star_indices);
    }

    indices
}

/// Stage 2:
/// This stage simply finds all neighbouring particles, using N_i(x*_i):
/// It returns a Vec, where Vec[i] is a vector of indices for neighbouring particles of particle[i]
/// This version uses a naive O(n^2) approach
#[inline(never)]
pub fn stage_2_naive(x_stars: &Vec<Vector3<f32>>) -> Vec<Vec<usize>> {
    //let mut indices = Vec::with_capacity(x_stars.len());
    let mut indices = Vec::new();

    for i in 0..x_stars.len() {
        let mut x_star_indices = Vec::new();
        for j in 0..x_stars.len() {
            if i != j && (x_stars[i] - x_stars[j]).norm_squared() <= RADIUS_H_POW2 {
                x_star_indices.push(j);
            }
        }
        indices.push(x_star_indices);
    }
    indices
}

/// Stage 3:
/// This stage performs all the updates which need to be done once per solver iteration:
/// 1) It calculates all of λ_i
/// 2) It calculates Δp_i, and then checks for continuous collision detection
/// 3) It updates all predicted positions: x*_i <= x*_i + Δp_i
/// NOTE: The way this function works is it mutates current predicted positions.
#[inline(never)]
fn stage_3(
    predicted_positions: &mut Vec<Vector3<f32>>,
    neighbour_indices: &Vec<Vec<usize>>,
    boundary: Boundary,
) {
    let mut lambdas = Vec::with_capacity(predicted_positions.len());

    // Calculate λ_i's
    for i in 0..predicted_positions.len() {
        let neighbours = &neighbour_indices[i];
        let lambda = lambda(&predicted_positions, neighbours, i);
        // println!("lambda {}", lambda);
        lambdas.push(lambda);
    }

    // Calculate Δp_i's
    for i in 0..predicted_positions.len() {
        let neighbours = &neighbour_indices[i];
        let delta_p = delta_p(&predicted_positions, neighbours, &lambdas, i);
        predicted_positions[i] += delta_p;

        predicted_positions[i].y = predicted_positions[i].y.clamp(0.0, f32::MAX);
        predicted_positions[i].x = predicted_positions[i].x.clamp(-boundary.x, boundary.x);
        predicted_positions[i].z = predicted_positions[i].z.clamp(-boundary.z, boundary.z);
    }
}

/// Stage 4:
/// This is the final stage of the algorithm.sat
/// 1) It updates velocity using a mixture of the old position and the predicted position:
/// v_i <= 1/Δt * (x*_i - x_i)
/// 2) It calculates (and applies) vorticity confinement and viscosity
/// 3) It assigns the final position: x_i <= x*_i
/// NOTE: This order slightly differs from the paper due to my implementation, however the end
/// result is the same.
#[inline(never)]
fn stage_4(
    particles: &mut Vec<Particle>,
    predicted_positions: Vec<Vector3<f32>>,
    neighbour_indices: Vec<Vec<usize>>,
    dt: f32,
) {
    let mut new_velocities = Vec::with_capacity(particles.len());

    // v_i <= 1/Δt * (x*_i - x_i)
    for i in 0..particles.len() {
        let new_velocity = (predicted_positions[i] - particles[i].position) / dt;
        new_velocities.push(new_velocity);
    }

    if std::env::args().find(|x| x == "nvort").is_none() {
        // First we apply vorticity confinement
        for i in 0..particles.len() {
            let vorticity = f_vorticity(
                &predicted_positions,
                &new_velocities,
                &neighbour_indices[i],
                i,
            );

            // F = ma, and m = 1 in this implementation
            // So F = a. Hence, we can directly update velocity like so:
            new_velocities[i] += dt * vorticity;
        }
    }

    if std::env::args().find(|x| x == "nvisc").is_none() {
        // Now we can apply viscosity
        for i in 0..particles.len() {
            new_velocities[i] = viscosity(
                &predicted_positions,
                &new_velocities,
                &neighbour_indices[i],
                i,
            );
        }
    }

    // Finally we can update position.
    // To do so, we'll just assign all newfound values into the original particles Vec
    for i in 0..particles.len() {
        particles[i].position = predicted_positions[i];
        particles[i].velocity = new_velocities[i];
    }
}

pub fn run_sim(
    particles: &mut Vec<Particle>,
    dt: f32,
    iterations: usize,
    boundary: Boundary,
    octree: bool,
) {
    let mut x_stars = stage_1(particles, dt);

    let neighbour_indices = if octree {
        stage_2_octree(&x_stars)
    } else {
        stage_2_naive(&x_stars)
    };

    for _ in 0..iterations {
        stage_3(&mut x_stars, &neighbour_indices, boundary);
    }

    stage_4(particles, x_stars, neighbour_indices, dt);

    if std::env::args().find(|x| x == "assert").is_some() {
        for particle in particles.iter() {
            assert!(!particle.position.x.is_nan());
            assert!(!particle.velocity.x.is_nan());
        }
    }
}
