mod single;
mod threaded;

pub use single::run_sim;
pub use threaded::run_sim_threaded;

// For benchmarking the octree
pub use single::{stage_2_octree, stage_2_naive};
pub use threaded::{stage_2_octree_threaded, stage_2_naive_threaded};
