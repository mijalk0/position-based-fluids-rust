use na::{IsometryMatrix3, Perspective3, Point3, Vector3};

pub struct Camera {
    pub position: Vector3<f32>,
    pub looking_at: Vector3<f32>,
    projection: [[f32; 4]; 4],
}

impl Camera {
    pub fn new(position: Vector3<f32>, looking_at: Vector3<f32>, aspect: f32) -> Self {
        Self {
            position,
            looking_at,
            projection: Perspective3::new(aspect, 0.7853982, 0.1, 300.0)
            //projection: Perspective3::new(aspect, 1.047198, 0.1, 100.0)
                .to_homogeneous()
                .into(),
        }
    }

    #[inline(always)]
    pub fn view_matrix(&self) -> [[f32; 4]; 4] {
        let eye = &Point3::from(self.position);
        let target = &Point3::from(self.looking_at);
        let up = &Vector3::y();

        IsometryMatrix3::look_at_rh(eye, target, up).to_homogeneous().into()
    }

    #[inline(always)]
    pub fn projection_matrix(&self) -> [[f32; 4]; 4] {
        self.projection
    }
}
