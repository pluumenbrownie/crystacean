from typing import Optional
import typer
from tqdm import tqdm
from typing_extensions import Annotated
import os
from shutil import copy2
from collections import Counter

import matplotlib.pyplot as plt
from matplotlib import rcParams
from ase.io import read as aseread

from classes import from_file as from_file_classes
from basis_vectors import full_lattice_from_basis_vectors
from crystacean import Lattice, from_dft_json
from cull_results import cull

# Test libraries required for other files.
import alive_progress, scipy


rcParams.update({"font.size": 11})
# lattice_size = 10
# site_size = 16
lattice_size = 6
site_size = 8
cross_thickness = 4

app = typer.Typer()


PLOT_BOOL = Annotated[
    bool,
    typer.Option(
        "--plot",
        "-p",
        help="Show the created base lattice and found interface configurations.",
    ),
]
DIRPATH_STR = Annotated[
    str,
    typer.Argument(
        help="The input directory to find interface configurations for.",
        show_default=False,
    ),
]
CDM_FLOAT = Annotated[
    float,
    typer.Option(
        "--creation-distance-margin",
        "-c",
        help="Maximum distance allowed between two attachment points when finding mid- and tripoints",
    ),
]
SAVE_DESCRIPTION = (
    "Save found interface configurations to given directory. Directory must not exist."
)
SAVETO_STR = Annotated[
    str,
    typer.Option("--save-to", "-s", help=SAVE_DESCRIPTION),
]
PARALLEL_BOOL = Annotated[
    bool,
    typer.Option(
        "--parallel", "-l", help="Use the parallel version of the solving method."
    ),
]
SINGLET_INT = Annotated[
    int,
    typer.Option(
        "--max-singlets",
        "-m",
        help="Maximum amount of singlets allowed for this new layer.",
    ),
]
SIMILARIRY_BOOL = Annotated[
    bool,
    typer.Option(
        "--similarity-filter",
        "-f",
        help="Use the similarity filter. This filter attempts to cull structures which are translations and rotations of existing structures.",
    ),
]
DD_FLOAT = Annotated[
    float,
    typer.Option(
        "--difference-distance",
        "-d",
        help="Maximum delta between distances allowed when using similarity filter. Does nothing when this filter is not used.",
    ),
]
RINGS_BOOL = Annotated[
    bool,
    typer.Option(
        "--rings-filter",
        "-r",
        help="Whether to apply the no-rings filter to the given lattice. Disable this when considering the first layer.",
    ),
]
ASE_BOOL = Annotated[
    bool,
    typer.Option(
        "--save-ase-json",
        "-j",
        help="If true, the found structures get exported as an ASE parsable structure in json format.",
    ),
]
STATISTICS = Annotated[
    bool,
    typer.Option(
        "--statistics",
        "-t",
        help="Prints how many found structures have the same amount of connections."
    )
]


@app.command()
def plot(path_to_json: Annotated[
    str,
    typer.Argument(help="Path to either a single json file, or a directory filled with json files.")
]):
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
    plot: PLOT_BOOL = False,
    save_to: SAVETO_STR = "",
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_BOOL = False,
    difference_distance: DD_FLOAT = 0.05,
):
    """
    Create SiO2 interface structures for a SiC unit cell with 4*size attachment points.
    """
    if not (plot or save_to):
        print("NOTE: both plot and save_to are false!")
    if save_to:
        os.mkdir(save_to)
    custom_boundary_points = full_lattice_from_basis_vectors(size)
    lattice = Lattice(custom_boundary_points)

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

    bit_lattice = lattice.get_intermediary(
        max_singlets=max_singlets,
        difference_distance=difference_distance,
        use_filter=similarity_filter,
    )
    print(bit_lattice)

    if use_parallel:
        solution = bit_lattice.solve_parallel(True)
    else:
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
            help="The input file to find interface configurations for. File must be in a ASE readable json format.",
            show_default=False,
        ),
    ],
    prefix: Annotated[
        str,
        typer.Argument(
            help="Identifying prefix to use when saving", show_default=False
        ),
    ],
    creation_distance_margin: CDM_FLOAT = 3.5,
    plot: PLOT_BOOL = False,
    save_to: SAVETO_STR = "",
    save_ase_json: ASE_BOOL = False,
    statistics: STATISTICS = False,
    rings_filter: RINGS_BOOL = False,
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_BOOL = False,
    difference_distance: DD_FLOAT = 0.05,
):
    """
    Create interface configurations from ASE json file. File must contain cell data.
    """
    if save_to:
        os.mkdir(save_to)
    ase_json_handler(
        filepath,
        prefix,
        creation_distance_margin,
        plot,
        save_to,
        save_ase_json,
        statistics,
        rings_filter,
        max_singlets,
        use_parallel,
        similarity_filter,
        difference_distance,
    )


@app.command()
def from_dft_folder(
    dirpath: DIRPATH_STR,
    save_to: Annotated[
        str,
        typer.Argument(help=SAVE_DESCRIPTION, show_default=False),
    ],
    prefix: Annotated[
        Optional[str], typer.Argument(help="Identifying prefix to use when saving", show_default="$dirpath lowest folder")
    ] = None,
    creation_distance_margin: CDM_FLOAT = 3.5,
    plot: PLOT_BOOL = False,
    statistics: STATISTICS = False,
    rings_filter: RINGS_BOOL = True,
    output_file_name: Annotated[
        str,
        typer.Option(help="Filename of cp2k DFT result file."),
    ] = "SiC-pos-1.xyz",
    test_mode: Annotated[
        bool,
        typer.Option("--test-mode", help="Stop after creating tempfile."),
    ] = False,
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_BOOL = False,
    difference_distance: DD_FLOAT = 0.05,
):
    """
    Find next layer configurations directly from CP2K DFT results.
    """
    assert os.path.isdir(dirpath), "Filepath should be a directory."
    dirpath = dirpath.removesuffix("/")
    save_to = save_to.removesuffix("/")
    if not prefix:
        prefix = dirpath.rsplit("/")[-1]

    try:
        os.mkdir(save_to)
    except FileExistsError:
        pass
    try:
        os.mkdir(save_to + "/temp")
    except FileExistsError:
        pass

    filepath = f"{save_to}/temp/{prefix}.json"

    cellfile = aseread(f"{dirpath}/new.xyz")
    cell = cellfile.get_cell()  # type: ignore
    tempfile = aseread(f"{dirpath}/{output_file_name}")
    tempfile.set_cell(cell)  # type: ignore
    tempfile.write(filepath)  # type: ignore

    if test_mode:
        return 0

    ase_json_handler(
        filepath,
        prefix,
        creation_distance_margin,
        plot,
        save_to,
        True,
        statistics,
        rings_filter,
        max_singlets,
        use_parallel,
        similarity_filter,
        difference_distance,
    )

    os.remove(f"{save_to}/temp/{prefix}.json")
    os.rmdir(f"{save_to}/temp")
    # if do_cp2kify:
    #     cp2kify(save_to, save_to)


@app.command()
def cp2kify(
    dirpath: DIRPATH_STR,
    save_to: Annotated[
        str,
        typer.Argument(help="Save found interface configurations to given directory.", show_default=False),
    ],
    in_path: Annotated[
        str,
        typer.Argument(
            help="The path to the `xxx.in` file needed for CP2K.", show_default=False
        ),
    ],
    basis_path: Annotated[
        str,
        typer.Argument(help="The `BASIS` file needed for CP2K.", show_default=False),
    ],
    run_path: Annotated[
        str,
        typer.Argument(
            help="The bash script to run the job with sbatch.", show_default=False
        ),
    ],
    destructive: Annotated[
        bool,
        typer.Option(
            "--destructive", help="Remove given json files when folder has been created."
        )
    ] = False,
):
    """
    Prepares a full folder of ASE json structures for batch processing with cp2k.
    """
    assert os.path.isdir(dirpath), "Filepath should be a directory."
    dirpath = dirpath.removesuffix("/")
    assert os.path.isfile(in_path), "xxx.in file not found."
    in_path = in_path.removesuffix("/")
    assert os.path.isfile(basis_path), "BASIS file not found."
    basis_path = basis_path.removesuffix("/")
    assert os.path.isfile(run_path), "Bash file not found."
    try:
         os.mkdir(save_to)
    except FileExistsError:
        pass
    save_to = save_to.removesuffix("/")

    file_list = tqdm(
        os.listdir(dirpath),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  ",
    )

    for structure_file in file_list:
        if not structure_file[-5:] == ".json":
            structure_file += ".json"
        os.mkdir(f"{save_to}/{structure_file[:-5]}")
        structure = aseread(f"{dirpath}/{structure_file}")
        structure.write(f"{save_to}/{structure_file[:-5]}/new.xyz")  # type: ignore
        copy2(in_path, f"{save_to}/{structure_file[:-5]}")
        copy2(basis_path, f"{save_to}/{structure_file[:-5]}")
        copy2(run_path, f"{save_to}/{structure_file[:-5]}")
        if destructive:
            os.remove(f"{dirpath}/{structure_file}")


@app.command()
def cull_results(
    dirpath: DIRPATH_STR,
    margin: Annotated[
        float,
        typer.Argument(
            help="Minimal deviation needed between structures to be designated as 'unique'.",
            show_default=False,
        ),
    ],
    postfix: Annotated[
        str,
        typer.Option(
            help="Identifier to indicate culled results. When empty, the margin is used."
        ),
    ] = "",
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
    creation_distance_margin: float,
    plot: bool,
    save_to: str,
    save_ase_json: bool,
    statistics: bool,
    rings_filter: bool,
    max_singlets: int,
    use_parallel: bool,
    similarity_filter: bool,
    difference_distance: float,
):
    if not (plot or save_to):
        print("NOTE: both plot and save_to are false!")
    lattice = from_dft_json(filepath, creation_distance_margin, False)

    if plot:
        plt.plot(*lattice.points_to_plot(), "o", markersize=lattice_size)
        plt.plot(
            *lattice.midpoints_to_plot(),
            "x",
            markersize=site_size,
            markeredgewidth=cross_thickness,
        )
        plt.plot(*lattice.tripoints_to_plot(), "s", markersize=site_size)
        plt.plot(*lattice.singlets_to_plot(), "^", markersize=site_size)
        plt.axis("equal")
        plt.xlabel("x (Å)")
        plt.ylabel("y (Å)")
        plt.show()
        # plt.plot(*lattice.points_to_plot(), "o", markersize=lattice_size)
        # plt.plot(
        #     *lattice.midpoints_to_plot(),
        #     "x",
        #     markersize=site_size,
        #     markeredgewidth=cross_thickness,
        # )
        # plt.plot(*lattice.tripoints_to_plot(), "s", markersize=site_size)
        # plt.axis("equal")
        # for num, x, y in lattice.no_rings_plot():
        #     plt.annotate(str(num), (x, y))
        # plt.show()

    bit_lattice = lattice.get_intermediary(
        max_singlets=max_singlets,
        difference_distance=difference_distance,
        use_filter=similarity_filter,
    )
    if rings_filter:
        noloops = lattice.no_rings()
        bit_lattice = bit_lattice.filtered(noloops)
    # print(bit_lattice)

    if use_parallel:
        solutions = bit_lattice.solve_parallel(True)
    else:
        solutions = bit_lattice.solve(True)

    progress = tqdm(
        enumerate(solutions),
        total=len(solutions),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  ",
    )

    if save_to:
        for number, solution in enumerate(solutions):
            solved_lattice = lattice.to_solved_lattice(solution)
            if save_ase_json:
                solved_lattice.export_as_ase_json(
                    f"{prefix}_{number:>04}.json", save_to
                )
            else:
                solved_lattice.export(save_to, f"{prefix}_{number:>04}.json")
    
    if statistics:
        solutions_dict = Counter()
        for solution in solutions:
            solved_lattice = lattice.to_solved_lattice(solution)
            type_amount = (
                len(solved_lattice.tripoints_to_plot()[0]),
                len(solved_lattice.midpoints_to_plot()[0]),
                len(solved_lattice.singlets_to_plot()[0]),
            )
            solutions_dict[type_amount] += 1
        for item in solutions_dict.items():
            print(item)

    if plot:
        for number, solution in progress:
            solved_lattice = lattice.to_solved_lattice(solution)
            progress.set_description(desc=f"Solution {number}")
            plt.plot(*solved_lattice.points_to_plot(), "o", markersize=lattice_size)
            #     for num, point in enumerate(zip( *lattice.points_to_plot())):
            #         plt.annotate(str(num), (point[0], point[1]))

            plt.plot(
                *solved_lattice.midpoints_to_plot(),
                "x",
                markersize=site_size,
                markeredgewidth=cross_thickness,
            )
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
