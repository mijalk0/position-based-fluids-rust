pub extern crate nalgebra as na;

pub mod artificial_pressure;
pub mod boundary;
pub mod camera;
pub mod incompressibility;
pub mod particle;
pub mod renderer;
pub mod simulation;
//pub mod spawn;
pub mod viscosity;
pub mod vorticity;

pub use boundary::Boundary;
pub use particle::Particle;
pub use simulation::{run_sim, run_sim_threaded};
//For benchmarking the octree
pub use simulation::{
    stage_2_naive, stage_2_naive_threaded, stage_2_octree, stage_2_octree_threaded,
};
pub use viscosity::viscosity;
pub use vorticity::f_vorticity;
