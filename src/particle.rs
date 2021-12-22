use na::Vector3;

/// Particles are really nothing more than position and velocity.
/// They all share mass and rest density. So there's not much setting them apart.
#[derive(Copy, Clone)]
pub struct Particle {
    pub position: Vector3<f32>,
    pub velocity: Vector3<f32>,
}

impl Particle {
    #[inline(always)]
    pub fn new(position: Vector3<f32>, velocity: Vector3<f32>) -> Self {
        Self { position, velocity }
    }
}
