#[macro_use]
extern crate criterion;

use criterion::{BatchSize, Criterion};
use final_report::*;
use nalgebra::Vector3;

fn generate_points(count: usize) -> Vec<Vector3<f32>> {
    let points: Vec<Vector3<f32>> = (0..count)
        .into_iter()
        .map(|_| {
            let mut vector = Vector3::from(rand::random::<[f32; 3]>());

            // Changes random samples from to bigger ranges
            vector.x *= 10.0;
            vector.y *= 5.0;
            vector.z *= 10.0;
            vector
        })
        .collect();

    points
}

fn generate_particles(count: usize) -> Vec<Particle> {
    let positions = generate_points(count);
    let velocities = vec![Vector3::from(rand::random::<[f32; 3]>()); count];

    let particles = positions
        .iter()
        .zip(velocities.iter())
        .map(|(position, velocity)| Particle::new(*position, *velocity))
        .collect();

    particles
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("simulation single-threaded octree", |b| {
        b.iter_batched(
            || generate_particles(5000),
            |mut particles| {
                let boundary = Boundary { x: 10.0, z: 10.0 };
                run_sim(&mut particles, 1.0 / 120.0, 3, boundary, true);
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("simulation multi-threaded octree", |b| {
        b.iter_batched(
            || generate_particles(5000),
            |mut particles| {
                let boundary = Boundary { x: 10.0, z: 10.0 };
                run_sim_threaded(&mut particles, 1.0 / 120.0, 3, boundary, true);
            },
            BatchSize::SmallInput,
        );
    });
    
    c.bench_function("simulation single-threaded naive", |b| {
        b.iter_batched(
            || generate_particles(5000),
            |mut particles| {
                let boundary = Boundary { x: 10.0, z: 10.0 };
                run_sim(&mut particles, 1.0 / 120.0, 3, boundary, false);
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("simulation multi-threaded naive", |b| {
        b.iter_batched(
            || generate_particles(5000),
            |mut particles| {
                let boundary = Boundary { x: 10.0, z: 10.0 };
                run_sim_threaded(&mut particles, 1.0 / 120.0, 3, boundary, false);
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("octree search single-threaded", |b| {
        b.iter_batched(
            || generate_points(5000),
            |points| {
                // Ignore return value
                let _ = stage_2_octree(&points);
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("octree search multi-threaded", |b| {
        b.iter_batched(
            || generate_points(5000),
            |points| {
                // Ignore return value
                let _ = stage_2_octree_threaded(&points);
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("naive search multi-threaded", |b| {
        b.iter_batched(
            || generate_points(5000),
            |points| {
                // Ignore return value
                let _ = stage_2_naive_threaded(&points);
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("naive search single-threaded", |b| {
        b.iter_batched(
            || generate_points(5000),
            |points| {
                // Ignore return value
                let _ = stage_2_naive(&points);
            },
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
