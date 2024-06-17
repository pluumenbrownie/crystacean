use ::lattice_solver::BitArrayFilter;
use ::lattice_solver::BitArraySettings;
use pyo3::prelude::*;
use std::ffi::OsString;

use ::lattice_solver::BitArrayRepresentation as WrappedRepresentation;
use ::lattice_solver::BitArraySolution as WrappedSolution;
use ::lattice_solver::Lattice as WrappedLattice;
use ::lattice_solver::site_filter::SiteFilter as WrappedFilter;

#[pyclass]
struct BitArraySolution {
    wrapped: WrappedSolution,
}

#[pymethods]
impl BitArraySolution {
    fn __str__(&self) -> String {
        self.wrapped.__str__()
    }
}

#[pyclass]
#[derive(Clone)]
struct SiteFilter {
    wrapped: WrappedFilter
}

impl SiteFilter {
    fn empty() -> Self {
        SiteFilter{wrapped: WrappedFilter::empty()}
    }
}

#[pyclass]
struct BitArrayRepresentation {
    wrapped: WrappedRepresentation,
}

#[pymethods]
impl BitArrayRepresentation {
    fn solve(&self, find_all: bool) -> Vec<BitArraySolution> {
        self.wrapped
            .solve(find_all, false)
            .into_iter()
            .map(|a| BitArraySolution { wrapped: a })
            .collect()
    }

    fn solve_parallel(&self, find_all: bool) -> Vec<BitArraySolution> {
        self.wrapped
            .solve_parallel(find_all, false)
            .into_iter()
            .map(|a| BitArraySolution { wrapped: a })
            .collect()
    }

    fn filtered(&self, filter: SiteFilter) -> BitArrayRepresentation {
        BitArrayRepresentation { wrapped: self.wrapped.filtered(filter.wrapped) }
    }

    fn __str__(&self) -> String {
        self.wrapped.__str__()
    }

    fn __repr__(&self) -> String {
        self.wrapped.__repr__()
    }
}

#[pyclass]
struct Lattice {
    wrapped: WrappedLattice,
}

#[pymethods]
impl Lattice {
    #[new]
    #[pyo3(signature = (input_lattice, distance_margin, autodetect_margin=true))]
    fn python_new(input_lattice: Vec<(Vec<f32>, Vec<Vec<f32>>)>, distance_margin: f32, autodetect_margin: bool) -> Self {
        Lattice {
            wrapped: WrappedLattice::python_new(input_lattice, distance_margin, autodetect_margin),
        }
    }

    fn points_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        self.wrapped.points_to_plot()
    }

    fn oxygens_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        self.wrapped.oxygens_to_plot()
    }

    fn tripoints_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        self.wrapped.tripoints_to_plot()
    }

    fn midpoints_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        self.wrapped.midpoints_to_plot()
    }

    fn singlets_to_plot(&self) -> (Vec<f32>, Vec<f32>) {
        self.wrapped.singlets_to_plot()
    }

    #[pyo3(signature = (max_singlets=2, difference_distance=0.05, max=None, use_filter=false))]
    fn get_intermediary(&self, max_singlets: usize, difference_distance: f32, max: Option<(f32, f32)>, use_filter: bool) -> BitArrayRepresentation {
        BitArrayRepresentation {
            wrapped: self.wrapped.get_intermediary(BitArraySettings::create(
                    max_singlets,
                    difference_distance,
                    max.unwrap_or(self.wrapped.find_max()),
                    if use_filter {
                        BitArrayFilter::Similarity
                    } else {
                        BitArrayFilter::None
                    }
                )
            ),
        }
    }

    fn no_rings(&self) -> SiteFilter {
        SiteFilter{wrapped: self.wrapped.no_rings()}
    }

    fn no_rings_plot(&self) -> Vec<(usize, f32, f32)>{
        self.wrapped.no_rings_plot()
    }

    fn to_solved_lattice(&self, solution: &BitArraySolution) -> Self {
        Self {
            wrapped: self.wrapped.to_solved_lattice(&solution.wrapped),
        }
    }

    pub fn export(&self, path: OsString, name: String) {
        self.wrapped.export(&path, name);
    }

    fn diagnostic_ase(&self) {
        self.wrapped.diagnostic_ase();
    }


    fn export_as_ase_json(&self, filename: String, folder: Option<String>) {
        if let Some(folder_path) = folder {
            let mut full_path = folder_path.clone();
            full_path.push('/');
            full_path.push_str(&filename);
            self.wrapped.export_as_ase_json(&full_path)
        } else {
            self.wrapped.export_as_ase_json(&filename)
        }
    }
}

#[pyfunction]
#[pyo3(signature = (filename, distance_margin, autodetect_margin=true))]
/// Put a new layer on top of an existing structure, calculated with DFT.
fn from_dft_json(filename: String, distance_margin: f32, autodetect_margin: bool) -> Lattice {
    Lattice {
        wrapped: WrappedLattice::from_dft_json(filename, distance_margin, autodetect_margin)
    }
}

#[pyfunction]
/// Test the import of the library
fn test_module() {
    println!("It works on my computer!")
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "crystacean")]
fn lattice_solver(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(test_module, m)?)?;
    m.add_function(wrap_pyfunction!(from_dft_json, m)?)?;
    m.add_class::<Lattice>()?;
    m.add_class::<BitArrayRepresentation>()?;
    m.add_class::<BitArraySolution>()?;
    Ok(())
}
