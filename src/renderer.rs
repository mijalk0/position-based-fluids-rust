use glium::index::PrimitiveType::TrianglesList;
use glium::{
    implement_vertex, uniform, Display, Frame, IndexBuffer, Program, Surface, VertexBuffer,
};
use obj::Obj;

use na::Vector3;

use crate::camera::Camera;
use crate::incompressibility::RHO_0;
use crate::particle::Particle;

const SCALING: f32 = 0.3 * RHO_0;

#[derive(Copy, Clone)]
struct Vertex {
    position: [f32; 3],
    normal: [f32; 3],
}

implement_vertex!(Vertex, position, normal);

#[derive(Copy, Clone)]
struct Movement {
    world_position: [f32; 3],
    velocity: [f32; 3],
    vel_color: [f32; 3],
}

implement_vertex!(Movement, world_position, velocity, vel_color);

pub struct SphereRenderer {
    shader: Program,
    vertices: VertexBuffer<Vertex>,
    indices: IndexBuffer<u16>,
    movement: VertexBuffer<Movement>,
}

impl SphereRenderer {
    pub fn new(display: &Display, obj: Obj, shader: Program, particles: &Vec<Particle>) -> Self {
        let mut vertices = Vec::new();

        for vertex in obj.vertices {
            let scaled_position = [
                vertex.position[0] * SCALING,
                vertex.position[1] * SCALING,
                vertex.position[2] * SCALING,
            ];

            let normal = na::Vector3::from(vertex.normal).normalize();
            vertices.push(Vertex {
                position: scaled_position,
                normal: <[f32; 3]>::from(normal),
            })
        }

        let mut indices = Vec::new();

        for index in obj.indices {
            indices.push(index);
        }

        let vertices = glium::VertexBuffer::new(display, &vertices).unwrap();
        let indices = glium::IndexBuffer::new(display, TrianglesList, &indices).unwrap();

        let mut movement = Vec::new();

        for _ in 0..particles.len() {
            movement.push(Movement {
                world_position: [0.0, 0.0, 0.0],
                velocity: [0.0, 0.0, 0.0],
                vel_color: [0.0, 0.0, 0.0],
            });
        }

        let movement = VertexBuffer::dynamic(display, &movement).unwrap();

        Self {
            shader,
            vertices,
            indices,
            movement,
        }
    }

    pub fn draw_particles(
        &mut self,
        surface: &mut Frame,
        particles: &Vec<Particle>,
        camera: &Camera,
    ) {
        for (particle, movement) in particles.iter().zip(self.movement.map().iter_mut()) {
            movement.world_position = <[f32; 3]>::from(particle.position);
            movement.velocity = <[f32; 3]>::from(particle.velocity);
            movement.vel_color = velocity_rgb(particle.velocity);
        }

        let parameters = glium::DrawParameters {
            backface_culling: glium::draw_parameters::BackfaceCullingMode::CullClockwise,
            depth: glium::Depth {
                test: glium::DepthTest::IfLessOrEqual,
                write: true,
                ..Default::default()
            },
            ..Default::default()
        };

        let uniforms = uniform! {
            view_matrix: camera.view_matrix(),
            proj_matrix: camera.projection_matrix(),
        };

        surface
            .draw(
                (&self.vertices, self.movement.per_instance().unwrap()),
                &self.indices,
                &self.shader,
                &uniforms,
                &parameters,
            )
            .unwrap();
    }
}

/// Converts velocity onto the rgb scale
#[inline(always)]
fn velocity_rgb(velocity: Vector3<f32>) -> [f32; 3] {
    const VELOCITY_THRESHHOLD: f32 = 20.0;
    const COLDEST_BLUE: f32 = 0.065;

    let speed =
        ((velocity.norm()).min(VELOCITY_THRESHHOLD) / VELOCITY_THRESHHOLD).max(COLDEST_BLUE);
    let color = colorgrad::turbo().at(speed.into()).to_linear_rgba();

    [color.0 as f32, color.1 as f32, color.2 as f32]
    // [0.8, 0.5, 0.0]
}
