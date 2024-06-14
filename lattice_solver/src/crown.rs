use std::f32::consts::PI;
use crate::*;

pub const SINGLE_CROWN: [f32; 3] = [0.0, 0.0, -1.7];

pub fn double_crown_rotated(theta: f32) -> [f32; 6] {
    // println!("{theta:?}");
    let r = (2.2252961107321143 - 1.0526814404964109_f32)
        .hypot(-1.3867293215799141 - -1.0917258490752173_f32);
    [
        r * (theta).cos(),
        r * (theta).sin(),
        -(9.26392293591872 - 8.41953003801284),
        r * (theta + PI).cos(),
        r * (theta + PI).sin(),
        -(9.26392293591872 - 8.313286837087986),
    ]
}

#[allow(clippy::suboptimal_flops)]
pub fn triple_crown_rotated(theta: f32) -> [f32; 9] {
    let r = (-0.000082766_f32).hypot(-1.397718937_f32);
    [
        r * (theta).cos(),
        r * (theta).sin(),
        -0.517599594,
        r * (theta + (2.0 / 3.0 * PI)).cos(),
        r * (theta + (2.0 / 3.0 * PI)).sin(),
        -0.517599594,
        r * (theta - (2.0 / 3.0 * PI)).cos(),
        r * (theta - (2.0 / 3.0 * PI)).sin(),
        -0.517599594,
    ]
}

pub fn double_angle(p1: &LatticePoint, p2: &LatticePoint) -> f32 {
    PI - ((p2.y - p1.y) / (p2.x - p1.x).hypot(p2.y - p1.y)).acos()
}
