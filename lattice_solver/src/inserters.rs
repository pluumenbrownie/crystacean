use std::collections::HashSet;
use kiddo::SquaredEuclidean;

use crate::*;

pub fn insert_singles(out_lattice: &mut Lattice) {
    for (number, silicon) in out_lattice
        .points
        .iter()
        .enumerate()
        .filter(|s| s.1.ghost_to.is_none())
    {
        out_lattice.oxygens.push(Oxygen::new(
            silicon.x,
            silicon.y,
            silicon.z - 1.7,
            SiteType::Singlet(Singlet([LatticeIndex(number)])),
        ));
    }
}

pub fn insert_midpoints(
    out_lattice: &mut Lattice,
    kdtree: &kiddo::float::kdtree::KdTree<f32, u64, 3, BINSIZE, u32>,
    node_search_distance: f32,
) {
    for (number, silicon) in out_lattice.points.iter().enumerate() {
        let mut close_points = kdtree
            .within::<SquaredEuclidean>(&[silicon.x, silicon.y, silicon.z], node_search_distance);
        // Sort results on lattice number
        close_points.sort_by_key(|p| p.item);
        let sites = close_points
            .iter()
            .skip(1)
            .filter(|s| s.item as usize > number)
            .filter(|s| {
                out_lattice.points[number].ghost_to.is_none()
                    || out_lattice.points[s.item as usize].ghost_to.is_none()
            });

        for site in sites {
            let x = (out_lattice.points[number].x + out_lattice.points[site.item as usize].x) / 2.0;
            let y = (out_lattice.points[number].y + out_lattice.points[site.item as usize].y) / 2.0;
            let z = (out_lattice.points[number].z + out_lattice.points[site.item as usize].z) / 2.0
                - 1.4;
            let sitetype = SiteType::Midpoint(Midpoint([
                LatticeIndex(number),
                LatticeIndex(site.item as usize),
            ]));
            out_lattice.oxygens.push(Oxygen::new(x, y, z, sitetype));
        }
    }
}

#[allow(clippy::float_cmp)]
pub fn insert_tripoints(
    out_lattice: &mut Lattice,
    kdtree: &kiddo::float::kdtree::KdTree<f32, u64, 3, BINSIZE, u32>,
    node_search_distance: f32,
) {
    let mut covered_sites = HashSet::new();
    for (number, silicon) in out_lattice.points.iter().enumerate() {
        let mut close_points = kdtree
            .within::<SquaredEuclidean>(&[silicon.x, silicon.y, silicon.z], node_search_distance);
        // Sort results on lattice number
        close_points.sort_by_key(|p| p.item);
        let sites = close_points
            .iter()
            .filter(|s| s.item as usize != number)
            .combinations(2)
            .filter(|a| {
                let mut identifier = [number as u64, a[0].item, a[1].item];
                identifier.sort_unstable();
                covered_sites.insert(identifier)
            })
            .filter(|a| {
                out_lattice.points[number].ghost_to.is_none()
                || out_lattice.points[a[0].item as usize].ghost_to.is_none()
                || out_lattice.points[a[1].item as usize].ghost_to.is_none()
            })
            .filter(|a| {
                out_lattice.points[a[0].item as usize]
                .distance_squared_to(&out_lattice.points[a[1].item as usize])
                <= node_search_distance
            })
            // .inspect(|a| if number == 25 {println!("({number}, {}, {})", a[0].item, a[1].item)})
            // .collect_vec()
            ;

        for site in sites {
            let x = (out_lattice.points[number].x
                + out_lattice.points[site[0].item as usize].x
                + out_lattice.points[site[1].item as usize].x)
                / 3.0;
            let y = (out_lattice.points[number].y
                + out_lattice.points[site[0].item as usize].y
                + out_lattice.points[site[1].item as usize].y)
                / 3.0;
            let z = (out_lattice.points[number].z
                + out_lattice.points[site[0].item as usize].z
                + out_lattice.points[site[1].item as usize].z)
                / 3.0
                - 1.1;
            let sitetype = SiteType::Tripoint(Tripoint([
                LatticeIndex(number),
                LatticeIndex(site[0].item as usize),
                LatticeIndex(site[1].item as usize),
            ]));
            out_lattice.oxygens.push(Oxygen::new(x, y, z, sitetype));
        }
    }
}
