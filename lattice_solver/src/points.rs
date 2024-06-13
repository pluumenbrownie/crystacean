use std::sync::{Arc, RwLock};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OxygenIndex(pub usize);

#[derive(Clone, Copy, Debug)]
pub struct LatticeIndex(pub usize);

#[derive(Debug)]
pub struct LatticePoint {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub connected_to: RwLock<Vec<OxygenIndex>>,
    pub ghost_to: Option<Arc<LatticePoint>>,
}

impl LatticePoint {
    pub fn new(x: f32, y: f32, z: f32, ghost_to: Option<Arc<Self>>) -> Arc<Self> {
        Arc::new(Self {
            x,
            y,
            z,
            connected_to: RwLock::new(vec![]),
            ghost_to,
        })
    }

    pub fn get_connections(&self) -> &RwLock<Vec<OxygenIndex>> {
        self.ghost_to
            .as_ref()
            .map_or(&self.connected_to, |point| &point.connected_to)
    }

    #[allow(clippy::suboptimal_flops)]
    pub fn distance_squared_to(&self, other: &Self) -> f32 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)
    }
}

#[derive(Clone)]
pub struct Oxygen {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub sitetype: SiteType,
    pub exclusions: Vec<OxygenIndex>,
}

impl Oxygen {
    pub const fn new(x: f32, y: f32, z: f32, sitetype: SiteType) -> Self {
        Self {
            x,
            y,
            z,
            sitetype,
            exclusions: vec![],
        }
    }
}

#[derive(Clone, Copy)]
pub enum SiteType {
    Singlet(Singlet),
    Midpoint(Midpoint),
    Tripoint(Tripoint),
}

impl SiteType {
    pub fn iter(&self) -> std::slice::Iter<'_, LatticeIndex> {
        match self {
            Self::Tripoint(c) => c.0.iter(),
            Self::Midpoint(c) => c.0.iter(),
            Self::Singlet(c) => c.0.iter(),
        }
    }
}

#[derive(Clone, Copy)]
pub struct Singlet(pub [LatticeIndex; 1]);

#[derive(Clone, Copy)]
pub struct Midpoint(pub [LatticeIndex; 2]);

#[derive(Clone, Copy)]
pub struct Tripoint(pub [LatticeIndex; 3]);
