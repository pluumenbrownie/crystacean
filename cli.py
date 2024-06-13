import typer
from tqdm import tqdm
from typing_extensions import Annotated
import os

import matplotlib.pyplot as plt
from matplotlib import rcParams
from ase.io import read as aseread

from classes import from_file as from_file_classes
from findthosepoints import full_lattice_from_basis_vectors
from crystacean import Lattice, from_dft_json
from cull_results import cull

# Test libraries required for other files.
import alive_progress, scipy


rcParams.update({"font.size": 11})
lattice_size = 10
site_size = 16
cross_thickness = 4

app = typer.Typer()


@app.command()
def plot(path_to_json: str):
    """
    Plot either a single or a directory of interface json files.
    """
    if os.path.isfile(path_to_json):
        try:
            lattice = from_file_classes(f"{path_to_json}")
            lattice.plot(draw_single=True)
        except KeyError:
            print(f"{path_to_json} has incorrect json layout for printing.")

    elif os.path.isdir(path_to_json):
        structures = [
            structure
            for structure in os.listdir(path_to_json)
            if structure[-5:] == ".json"
        ]
        print(f"Found {len(structures)} json files in directory.")
        for structure in structures:
            try:
                lattice = from_file_classes(f"{path_to_json}/{structure}")
                print(structure)
                lattice.plot(draw_single=True)
            except KeyError:
                print(f"{structure} has incorrect json layout for printing.")

    else:
        raise FileNotFoundError(f"File or directory {path_to_json} not found.")


@app.command()
def from_size(
    name: str,
    size: Annotated[
        int,
        typer.Option(
            help="Size of base lattice to use. Amount of attachment sites is 4*size."
        ),
    ] = 1,
    plot: Annotated[
        bool,
        typer.Option(
            help="Show the created base lattice and found interface configurations."
        ),
    ] = False,
    save_to: Annotated[
        str,
        typer.Option(
            help="Save found interface configurations to given directory. Directory must not exist."
        ),
    ] = "",
):
    """
    Create SiO2 interface structures for a SiC unit cell with 4*size attachment points.
    """
    if not (plot or save_to):
        print("NOTE: both plot and save_to are false!")
    if save_to:
        os.mkdir(save_to)
    custom_boundary_points = full_lattice_from_basis_vectors(size)
    lattice = Lattice(custom_boundary_points, 1.1)

    if plot:
        plt.plot(
            *lattice.points_to_plot(),
            "o",
            markersize=lattice_size,
            label="Previous layer",
        )
        plt.plot(
            *lattice.midpoints_to_plot(),
            "x",
            markersize=site_size,
            label="Midpoint",
            markeredgewidth=4,
        )
        plt.plot(
            *lattice.tripoints_to_plot(), "s", markersize=site_size, label="Tripoint"
        )
        plt.plot(*lattice.singlets_to_plot(), "^", markersize=site_size, label="Single")
        plt.axis("equal")
        plt.xlabel("x (Å)")
        plt.ylabel("y (Å)")
        plt.legend(
            loc="upper center", bbox_to_anchor=(0.5, 1.11), ncol=4, fancybox=True
        )
        plt.show()

    bit_lattice = lattice.get_intermediary()
    print(bit_lattice)

    solutions = bit_lattice.solve(True)

    progress = tqdm(
        enumerate(solutions),
        total=len(solutions),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  ",
    )
    for number, solution in progress:
        solved_lattice = lattice.to_solved_lattice(solution)
        if save_to:
            solved_lattice.export(save_to, f"{name}_{number:>06}.json")

        if plot:
            progress.set_description(desc=f"Solution {number}")
            plt.plot(
                *solved_lattice.points_to_plot(),
                "o",
                markersize=lattice_size,
                label="Previous layer",
            )
            #     for num, point in enumerate(zip( *lattice.points_to_plot())):
            #         plt.annotate(str(num), (point[0], point[1]))

            plt.plot(
                *solved_lattice.midpoints_to_plot(),
                "x",
                markersize=site_size,
                label="Midpoint", 
                markeredgewidth=cross_thickness,
            )
            plt.plot(
                *solved_lattice.tripoints_to_plot(),
                "s",
                markersize=site_size,
                label="Tripoint",
            )
            plt.plot(
                *solved_lattice.singlets_to_plot(),
                "^",
                markersize=site_size,
                label="Single",
            )
            #     for num, point in enumerate(zip( *solved_lattice.oxygens_to_plot())):
            #         plt.annotate(str(num), (point[0], point[1]))
            plt.axis("equal")
            plt.xlabel("x (Å)")
            plt.ylabel("y (Å)")
            plt.legend(
                loc="upper center", bbox_to_anchor=(0.5, 1.11), ncol=4, fancybox=True
            )
            plt.show()


@app.command()
def from_file(
    filepath: Annotated[
        str,
        typer.Argument(
            help="The input file to find interface configurations for."
        )
    ],
    prefix: Annotated[
        str,
        typer.Argument(
            help="Identifying prefix to use when saving"
        )
    ],
    distance_margin: Annotated[
        float,
        typer.Option(
            help="Maximum distance allowed between two attachment points when finding mid- and tripoints"
        ),
    ] = 3.5,
    plot: Annotated[
        bool,
        typer.Option(
            help="Show the created base lattice and found interface configurations."
        ),
    ] = False,
    save_to: Annotated[
        str,
        typer.Option(
            help="Save found interface configurations to given directory. Directory must not exist."
        ),
    ] = "",
    filtered: Annotated[
        bool,
        typer.Option(
            help="Whether to apply the no-rings filter to the given lattice. Disable this when considering the first layer."
        ),
    ] = True,
):
    """
    Create interface configurations from ASE json file. File must contain cell data.
    """
    if save_to:
        os.mkdir(save_to)
    ase_json_handler(filepath, prefix, distance_margin, plot, save_to, filtered)


@app.command()
def from_dft_folder(
    dirpath: Annotated[
        str,
        typer.Argument(
            help="The input directory to find interface configurations for."
        )
    ],
    save_to: Annotated[
        str,
        typer.Argument(
            help="Save found interface configurations to given directory. Directory must not exist."
        ),
    ],
    prefix: Annotated[
        str,
        typer.Argument(
            help="Identifying prefix to use when saving"
        )
    ],
    distance_margin: Annotated[
        float,
        typer.Option(
            help="Maximum distance allowed between two attachment points when finding mid- and tripoints"
        ),
    ] = 3.5,
    plot: Annotated[
        bool,
        typer.Option(
            help="Show the created base lattice and found interface configurations."
        ),
    ] = False,
    filtered: Annotated[
        bool,
        typer.Option(
            help="Whether to apply the no-rings filter to the given lattice. Disable this when considering the first layer."
        ),
    ] = True,
    output_file_name: Annotated[
        str,
        typer.Option(
            help="Filename of cp2k DFT result file."
        ),
    ] = "SiC-pos-1.xyz",
    test_mode: Annotated[
        bool,
        typer.Option(
            help="Stop after creating tempfile."
        ),
    ] = False
):
    """
    Find next layer configurations directly from CP2K DFT results.
    """
    assert os.path.isdir(dirpath), "Filepath should be a directory."
    dirpath.removesuffix("/")
    save_to.removesuffix("/")

    try:
        os.mkdir(save_to)
    except FileExistsError:
        pass
    try:
        os.mkdir(save_to + "/temp")
    except FileExistsError:
        pass

    cellfile = aseread(f"{dirpath}/new.xyz")
    cell = cellfile.get_cell() # type: ignore
    tempfile = aseread(f"{dirpath}/{output_file_name}")
    tempfile.set_cell(cell) # type: ignore
    tempfile.write(f"{save_to}/temp/{prefix}.json") # type: ignore

    if test_mode:
        return 0
    
    filepath = f"{save_to}/temp/{prefix}.json"
    ase_json_handler(filepath, prefix, distance_margin, plot, save_to, filtered)

    os.remove(f"{save_to}/temp/{prefix}.json")
    os.rmdir(f"{save_to}/temp")


@app.command()
def cull_results(
    dirpath: Annotated[
        str,
        typer.Argument(
            help="The input directory to find interface configurations for."
        )
    ],
    margin: Annotated[
        float,
        typer.Argument(
            help="Minimal deviation needed between structures to be designated as 'unique'."
        ),
    ],
    postfix: Annotated[
        str,
        typer.Option(
            help="Identifier to indicate culled results. When empty, the margin is used."
        )
    ] = ""
):
    """
    Remove structures which are translated and/or rotated duplicates of existing ones.
    """
    if not postfix:
        postfix = str(margin).replace(".", "_")
    cull(dirpath, margin, postfix)

    

def ase_json_handler(
    filepath: str,
    prefix: str,
    distance_margin: float,
    plot: bool,
    save_to: str,
    filtered: bool,
):
    if not (plot or save_to):
        print("NOTE: both plot and save_to are false!")
    lattice = from_dft_json(filepath, distance_margin, False)

    if plot:
        plt.plot(*lattice.points_to_plot(), "o", markersize=lattice_size)
        plt.plot(*lattice.midpoints_to_plot(), "x", markersize=site_size, markeredgewidth=cross_thickness)
        plt.plot(*lattice.tripoints_to_plot(), "s", markersize=site_size)
        plt.plot(*lattice.singlets_to_plot(), "^", markersize=site_size)
        plt.axis("equal")
        plt.xlabel("x (Å)")
        plt.ylabel("y (Å)")
        plt.show()

    bit_lattice = lattice.get_intermediary()
    if filtered:
        noloops = lattice.no_rings()
        bit_lattice = bit_lattice.filtered(noloops)
    print(bit_lattice)

    solutions = bit_lattice.solve(True)

    progress = tqdm(
        enumerate(solutions),
        total=len(solutions),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  ",
    )

    for number, solution in progress:
        solved_lattice = lattice.to_solved_lattice(solution)
        if save_to:
            solved_lattice.export(save_to, f"{prefix}_{number:>04}.json")

        if plot:
            progress.set_description(desc=f"Solution {number}")
            plt.plot(*solved_lattice.points_to_plot(), "o", markersize=lattice_size)
            #     for num, point in enumerate(zip( *lattice.points_to_plot())):
            #         plt.annotate(str(num), (point[0], point[1]))

            plt.plot(*solved_lattice.midpoints_to_plot(), "x", markersize=site_size, markeredgewidth=cross_thickness)
            plt.plot(*solved_lattice.tripoints_to_plot(), "s", markersize=site_size)
            plt.plot(*solved_lattice.singlets_to_plot(), "^", markersize=site_size)
            #     for num, point in enumerate(zip( *solved_lattice.oxygens_to_plot())):
            #         plt.annotate(str(num), (point[0], point[1]))
            plt.axis("equal")
            plt.xlabel("x (Å)")
            plt.ylabel("y (Å)")
            plt.show()


if __name__ == "__main__":
    app()
