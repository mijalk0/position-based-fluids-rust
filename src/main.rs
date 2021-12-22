extern crate final_report;
extern crate glium;
extern crate rand;
extern crate rand_chacha;

use crate::glutin::dpi::PhysicalPosition;
use final_report::{
    camera::Camera, particle::Particle, renderer::SphereRenderer, run_sim, run_sim_threaded, Boundary,
};
use glium::{glutin, Program, Surface};
use nalgebra::Vector3;
use obj::{load_obj, Obj};
use rand::{Rng, SeedableRng};
use sphrs::{Coordinates, SHCoordinates};
use std::f32::consts::PI;
use std::fs::File;
use std::io::BufReader;
use std::sync::mpsc::channel;
use std::thread;
use std::time::{Duration, Instant};

fn main() {
    // Boilerplate glium setup
    let event_loop = glutin::event_loop::EventLoop::new();
    let wb = glutin::window::WindowBuilder::new().with_title("Position Based Fluids");
    let cb = glutin::ContextBuilder::new();
    let display = glium::Display::new(wb, cb, &event_loop).unwrap();

    // Initialize random number generator for y values
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(0);

    let mut particles = Vec::new();

    let n = extract_usize("n", 3000);

    let x_boundary = extract_f32("x", 10.0);
    let y_range = extract_f32("y", 10.0);
    let z_boundary = extract_f32("z", 10.0);

    // Totally random locations
    for _ in 0..n {
        let pos_x = rng.gen_range(-x_boundary..x_boundary);
        let pos_y = rng.gen_range(0.0..y_range);
        let pos_z = rng.gen_range(-z_boundary..z_boundary);

        let position = Vector3::new(pos_x, pos_y, pos_z);
        let velocity = Vector3::zeros();

        particles.push(Particle::new(position, velocity));
    }

    let camera_position = Vector3::new(20.0, 30.0, 40.0);
    let looking_at = Vector3::zeros();

    let width = display.gl_window().window().inner_size().width as f32;
    let height = display.gl_window().window().inner_size().height as f32;
    let aspect = width / height;

    let mut camera = Camera::new(camera_position, looking_at, aspect);

    let sphere_path = if std::env::args().find(|x| x == "smooth").is_some() {
        "sphere_smooth.obj"
    } else {
        "sphere.obj"
    };

    let octree = std::env::args().find(|x| x == "noctree").is_none();

    // Need to load the sphere
    let sphere_buf = BufReader::new(File::open(sphere_path).unwrap());
    let sphere_obj: Obj = load_obj(sphere_buf).unwrap();

    // Load up our shaders
    let vertex_src = include_str!("sphere.vs");
    let frag_src = include_str!("sphere.fs");
    let shader = Program::from_source(&display, vertex_src, frag_src, None).unwrap();

    // Create our sphere renderer
    let mut sphere_renderer = SphereRenderer::new(&display, sphere_obj, shader, &particles);

    // Setting up the physics to run on another thread
    let (send, recv) = channel();

    let particles_cloned = particles.clone();
    let boundary = Boundary {
        x: x_boundary,
        z: z_boundary,
    };

    thread::spawn(move || {
        let mut particles = particles_cloned;

        // Aiming for 60fps (2 substep) animation ...?
        let mut dt = 1.0 / 60.0 / 2.0;
        let mut time_counter = Instant::now();
        std::thread::sleep(Duration::from_millis(100));

        // let mut iteration = 0;
        loop {
            if std::env::args().find(|x| x == "realtime").is_some() {
                dt = (Instant::now() - time_counter).as_secs_f32();
                time_counter = Instant::now();
            }

            // Handle simulation
            //println!("iteration {}", iteration);
            
            if std::env::args().find(|x| x == "nthreaded").is_some() {
                run_sim(&mut particles, dt, 3, boundary, octree);
            } else {
                run_sim_threaded(&mut particles, dt, 3, boundary, octree);
            }
            //println!("particle position: {:#}", particles[0].position);
            //println!("particle velocity: {:#}", particles[0].velocity);
            // iteration += 1;
            // println!("iteration {}", iteration);
            send.send(particles.clone()).unwrap();
        }
    });

    print_usage();

    let mut left_mouse_pressed = false;
    let mut curr_pixels: Option<glutin::dpi::PhysicalPosition<f64>> = None;

    event_loop.run(move |ev, _, control_flow| {
        let next_frame_time =
            std::time::Instant::now() + std::time::Duration::from_nanos(16_666_667);

        // This is just glium boiler plate to handle window close events and frame time
        *control_flow = glutin::event_loop::ControlFlow::WaitUntil(next_frame_time);
        match ev {
            glutin::event::Event::WindowEvent { event, .. } => match event {
                glutin::event::WindowEvent::CloseRequested => {
                    *control_flow = glutin::event_loop::ControlFlow::Exit;
                    return;
                }
                glutin::event::WindowEvent::MouseInput { button, state, .. } => match button {
                    glutin::event::MouseButton::Left => match state {
                        glutin::event::ElementState::Pressed => {
                            left_mouse_pressed = true;
                        }
                        glutin::event::ElementState::Released => {
                            left_mouse_pressed = false;
                            curr_pixels = None;
                        }
                    },
                    _ => return,
                },
                glutin::event::WindowEvent::CursorMoved { position, .. } => {
                    if left_mouse_pressed {
                        camera_orbit(&mut camera, &mut curr_pixels, position);
                    }
                }
                _ => return,
            },
            glutin::event::Event::DeviceEvent { event, .. } => match event {
                glutin::event::DeviceEvent::MouseWheel { delta } => match delta {
                    glutin::event::MouseScrollDelta::PixelDelta(pixel_delta) => {
                        camera_zoom(&mut camera, pixel_delta);
                    }
                    _ => return,
                },
                _ => return,
            },
            _ => (),
        }

        if let Ok(new_particles) = recv.try_recv() {
            particles = new_particles;
        }

        // Drawing
        let mut target = display.draw();
        target.clear_color_and_depth((1.0, 1.0, 0.8, 1.0), 1.0);
        sphere_renderer.draw_particles(&mut target, &particles, &camera);
        target.finish().unwrap();
    });
}

fn print_usage() {
    println!("---------------------------------Viewing---------------------------------");
    println!("Scroll up and down to zoom in and out");
    println!("Click and drag to orbit around the origin (0, 0, 0)");
    println!("------------------------------Possible Args------------------------------");
    println!("\"smooth\" to run with smoother sphere meshes");
    println!("\"realtime\" to run with realtime time skips (unstable), default time step is 1.0/120.0 seconds per update");
    println!("\"assert\" to check for NaN occurences");
    println!("\"nthreaded\" to run without multithreading implementation");
    println!("\"npress\" to run with no artificial pressure term");
    println!("\"nvisc\" to run with no XSPH viscosity");
    println!("\"nvort\" to run with no vorticity confinement");
    println!("\"noctree\" to run neighbour search without octree (use O(n^2) method instead)");
    println!("---------------------------Optional Parameters---------------------------");
    println!("\"n #\" to run with # particles (default 3000)");
    println!("\"y #\" to spawn particles within [0,#] (default 10.0)");
    println!("\"x #\" to run with an x boundary of [-#,#] (default 10.0)");
    println!("\"z #\" to run with an z boundary of [-#,#] (default 10.0)");
}

fn camera_zoom(camera: &mut Camera, pixel_delta: glutin::dpi::PhysicalPosition<f64>) {
    const ZOOM_SCALE: f32 = 0.02;
    const CAMERA_MIN_ZOOM: f32 = 15.0;
    const CAMERA_MAX_ZOOM: f32 = 150.0;
    let scaled_delta =
        (pixel_delta.y as f32 * ZOOM_SCALE) * (camera.position - camera.looking_at).normalize();
    if CAMERA_MIN_ZOOM <= camera.position.norm() && camera.position.norm() <= CAMERA_MAX_ZOOM {
        let new_position = camera.position + scaled_delta;
        // Clamp within range (A litt more tight for float safety)

        if new_position.normalize().x.signum() == camera.position.normalize().x.signum() {
            camera.position = new_position;
            let camera_clamped_magnitude = camera
                .position
                .norm()
                .clamp(CAMERA_MIN_ZOOM + 1.0, CAMERA_MAX_ZOOM - 1.0);

            camera.position = camera.position.normalize() * camera_clamped_magnitude;
        }
    }
}

fn camera_orbit(
    camera: &mut Camera,
    curr_pixels: &mut Option<PhysicalPosition<f64>>,
    position: PhysicalPosition<f64>,
) {
    const Y_SCALE: f32 = 0.005;
    const X_SCALE: f32 = 0.005;

    let curr_position: Coordinates<f32> =
        Coordinates::cartesian(camera.position.z, camera.position.x, camera.position.y);

    let (theta, phi) = if let Some(pixels) = curr_pixels {
        let delta = (pixels.y - position.y) as f32 * Y_SCALE;
        let mut theta = curr_position.theta() + delta;
        theta = theta.clamp(PI / 6.0, PI - PI / 6.0);

        let delta = (pixels.x - position.x) as f32 * X_SCALE;
        let phi = curr_position.phi() + delta;

        *curr_pixels = Some(position);
        (theta, phi)
    } else {
        *curr_pixels = Some(position);
        (curr_position.theta(), curr_position.phi())
    };

    let new_position = Coordinates::spherical(curr_position.r(), theta, phi);
    camera.position = Vector3::new(new_position.y(), new_position.z(), new_position.x());
}

fn extract_f32(key: &str, default: f32) -> f32 {
    if let Some(i) = std::env::args().position(|x| x == key) {
        std::env::args()
            .collect::<Vec<String>>()
            .get(i + 1)
            .expect(format!("No number specified for {}", key).as_str())
            .parse::<f32>()
            .expect(format!("No number specified for {}", key).as_str())
    } else {
        default
    }
}

fn extract_usize(key: &str, default: usize) -> usize {
    if let Some(i) = std::env::args().position(|x| x == key) {
        std::env::args()
            .collect::<Vec<String>>()
            .get(i + 1)
            .expect(format!("No number specified for {}", key).as_str())
            .parse::<usize>()
            .expect(format!("No number specified for {}", key).as_str())
    } else {
        default
    }
}
