from typing import Optional
import typer
from tqdm import tqdm
from typing_extensions import Annotated
import os
from shutil import copy2
from collections import Counter
import json
import re

import matplotlib.pyplot as plt
from matplotlib import rcParams
from ase.io import read as aseread

from classes import from_file as from_file_classes
from basis_vectors import full_lattice_from_basis_vectors
from crystacean import Lattice, from_dft_json  # type: ignore
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


def is_file_callback(path: str):
    if not os.path.isfile(path):
        raise typer.BadParameter(f"{path} is not a file.")
    return path


def is_dir_callback(path: str):
    if not os.path.isdir(path):
        raise typer.BadParameter(f"{path} is not a directory.")
    path.removesuffix("/")
    return path


def path_exists_callback(path: str):
    if not os.path.exists(path):
        raise typer.BadParameter(f"{path} does not exist.")
    return path


def path_empty_callback(path: str):
    if os.path.exists(path):
        raise typer.BadParameter(f"{path} already exists.")
    return path


HELP_CREATE = "Creation"
HELP_INPUT = "Input"
HELP_OUTPUT = "Output"
HELP_DEBUG = "Debugging"

PLOT_BOOL = Annotated[
    bool,
    typer.Option(
        "--plot",
        "-p",
        help="Show the created base lattice and found interface configurations.",
        rich_help_panel=HELP_OUTPUT,
    ),
]
DIRPATH_STR = Annotated[
    str,
    typer.Argument(
        help="The input directory to find interface configurations for.",
        show_default=False,
        callback=is_dir_callback,
    ),
]
CDM_FLOAT = Annotated[
    float,
    typer.Option(
        "--creation-distance-margin",
        "-c",
        help="Maximum distance allowed between two attachment points when finding mid- and tripoints",
        rich_help_panel=HELP_CREATE,
    ),
]
SAVE_DESCRIPTION = (
    "Save found interface configurations to given directory. Directory must not exist."
)
SAVETO_STR = Annotated[
    str,
    typer.Option(
        "--save-to",
        "-s",
        help=SAVE_DESCRIPTION,
        callback=path_empty_callback,
        rich_help_panel=HELP_OUTPUT,
    ),
]
PARALLEL_BOOL = Annotated[
    bool,
    typer.Option(
        "--parallel",
        "-l",
        help="Use the parallel version of the solving method. Disabled due to lack of utility.",
        hidden=True,
        rich_help_panel=HELP_CREATE,
    ),
]
SINGLET_INT = Annotated[
    int,
    typer.Option(
        "--max-singlets",
        "-m",
        help="Maximum amount of singlets allowed for this new layer.",
        rich_help_panel=HELP_CREATE,
    ),
]
SIMILARIRY_OPT_FLOAT = Annotated[
    Optional[float],
    typer.Option(
        "--similarity-filter",
        "-f",
        help="Use the similarity filter with the given margin. This filter attempts to cull structures which are translations and rotations of existing structures.",
        rich_help_panel=HELP_CREATE,
        show_default=False,
    ),
]
RINGS_BOOL = Annotated[
    bool,
    typer.Option(
        "--rings-filter",
        "-r",
        help="Apply the no-rings filter to the given lattice. Enable this when considering layers beyond the first.",
        rich_help_panel=HELP_CREATE,
    ),
]
ASE_BOOL = Annotated[
    bool,
    typer.Option(
        "--save-ase-json",
        "-j",
        help="If true, the found structures get exported as an ASE parsable structure in json format. Only works with save-to (-s) set.",
        rich_help_panel=HELP_OUTPUT,
    ),
]
STATISTICS_BOOL = Annotated[
    bool,
    typer.Option(
        "--statistics",
        "-t",
        help="Prints how many found structures have the same amount of connections.",
        rich_help_panel=HELP_OUTPUT,
    ),
]


@app.command()
def plot(
    path_to_json: Annotated[
        str,
        typer.Argument(
            help="Path to either a single json file, or a directory filled with json files.",
            callback=path_exists_callback,
        ),
    ]
):
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
            help="Size of base lattice to use. Amount of attachment sites is 4*size.",
            rich_help_panel=HELP_CREATE,
        ),
    ] = 1,
    plot: PLOT_BOOL = False,
    save_to: SAVETO_STR = "",
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_OPT_FLOAT = None,
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

    if similarity_filter:
        difference_distance = similarity_filter
    else:
        difference_distance = 0.0
    bit_lattice = lattice.get_intermediary(
        max_singlets=max_singlets,
        difference_distance=difference_distance,
        use_filter=similarity_filter is not None,
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
    rings_filter: RINGS_BOOL = False,
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_OPT_FLOAT = None,
    plot: PLOT_BOOL = False,
    statistics: STATISTICS_BOOL = False,
    save_to: SAVETO_STR = "",
    save_ase_json: ASE_BOOL = False,
    debug: Annotated[
        bool,
        typer.Option(
            "--debug", help="Plot debug information.", rich_help_panel=HELP_DEBUG
        ),
    ] = False,
    silent: Annotated[
        bool,
        typer.Option(
            "--silent", help="Don't show the progress bar.", rich_help_panel=HELP_OUTPUT
        ),
    ] = False,
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
        debug,
        silent,
    )


@app.command()
def from_dft_folder(
    dirpath: DIRPATH_STR,
    save_to: Annotated[
        str,
        typer.Argument(help=SAVE_DESCRIPTION, show_default=False),
    ],
    prefix: Annotated[
        Optional[str],
        typer.Argument(
            help="Identifying prefix to use when saving",
            show_default="$dirpath lowest folder",
        ),
    ] = None,
    creation_distance_margin: CDM_FLOAT = 3.5,
    max_singlets: SINGLET_INT = 2,
    use_parallel: PARALLEL_BOOL = False,
    similarity_filter: SIMILARIRY_OPT_FLOAT = None,
    statistics: STATISTICS_BOOL = False,
    plot: PLOT_BOOL = False,
    output_file_name: Annotated[
        str,
        typer.Option(
            help="Filename of cp2k DFT result file.", rich_help_panel=HELP_INPUT
        ),
    ] = "SiC-pos-1.xyz",
    test_mode: Annotated[
        bool,
        typer.Option(
            "--test-mode",
            help="Stop after creating tempfile.",
            rich_help_panel=HELP_DEBUG,
        ),
    ] = False,
    disable_rings_filter: Annotated[
        bool,
        typer.Option(
            "--disable-rings-filter",
            help="Disable the rings filter when generating new structures.",
            rich_help_panel=HELP_CREATE,
        ),
    ] = False,
):
    """
    Find next layer configurations directly from CP2K DFT results.
    """
    save_to = save_to.removesuffix("/")
    if not prefix:
        dirpath = dirpath.strip("/")
        prefix = dirpath.split("/")[-1]

    try:
        os.mkdir(save_to)
    except FileExistsError:
        pass
    try:
        os.mkdir(save_to + "/temp")
    except FileExistsError:
        pass

    filepath = f"{save_to}/temp/{prefix}.json"

    cellfile = aseread(f"{dirpath}/SiC-1.restart")
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
        not disable_rings_filter,
        max_singlets,
        use_parallel,
        similarity_filter,
        False,
        False,
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
        typer.Argument(
            help="Save found interface configurations to given directory.",
            show_default=False,
        ),
    ],
    in_path: Annotated[
        str,
        typer.Argument(
            help="The path to the `xxx.in` file needed for CP2K.",
            show_default=False,
            callback=is_file_callback,
        ),
    ],
    basis_path: Annotated[
        str,
        typer.Argument(
            help="The `BASIS` file needed for CP2K.",
            show_default=False,
            callback=is_file_callback,
        ),
    ],
    run_path: Annotated[
        str,
        typer.Argument(
            help="The bash script to run the job with sbatch.",
            show_default=False,
            callback=is_file_callback,
        ),
    ],
    destructive: Annotated[
        bool,
        typer.Option(
            "--destructive",
            help="Remove given json files when folder has been created.",
        ),
    ] = False,
):
    """
    Prepares a full folder of ASE json structures for batch processing with cp2k.
    """
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


@app.command()
def process_results(
    dirpath: DIRPATH_STR,
    type_info_path: Annotated[
        str,
        typer.Argument(
            help="The directory with the original ASE json files, which contain the amount "
            + "of tripoints, midpoints and singlets for a given structure.",
            callback=is_dir_callback,
            show_default=False,
        ),
    ],
    output_file_path: Annotated[
        str,
        typer.Argument(help="File to store the found results in.", show_default=False),
    ],
):
    """
    Extract energies from DFT results and write them to a file.
    """
    folder_list = tqdm(
        os.listdir(dirpath),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}  ",
    )

    if not output_file_path[-4:] == ".csv":
        output_file_path += ".csv"
    with open(output_file_path, "x") as output_file:
        output_file.write("structure,energy,tripoints,midpoints,singlets\n")
        for folder in folder_list:
            if os.path.isfile(folder):
                continue
            structure_id = "-".join(re.findall(r"((?:\d{4,})+)(?!.out)", folder))
            if structure_id is None:
                continue
            structure_file = aseread(f"{dirpath}/{folder}/SiC-pos-1.xyz")
            energy = structure_file.info["E"]  # type: ignore

            with open(f"{type_info_path}/{folder}.json") as file:
                type_info = json.load(file)
            tripoints = type_info["1"]["tripoints"]
            midpoints = type_info["1"]["midpoints"]
            singlets = type_info["1"]["singlets"]

            output_file.write(
                f"{structure_id},{energy},{tripoints},{midpoints},{singlets}\n"
            )


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
    similarity_filter: Optional[float],
    debug: bool,
    silent: bool,
):
    if not (plot or save_to):
        print("NOTE: both plot and save_to are false!")
    lattice = from_dft_json(filepath, creation_distance_margin, False)

    if debug:
        plt.plot(*lattice.points_to_plot(), "o", markersize=lattice_size)
        plt.plot(
            *lattice.midpoints_to_plot(),
            "x",
            markersize=site_size,
            markeredgewidth=cross_thickness,
        )
        plt.plot(*lattice.tripoints_to_plot(), "s", markersize=site_size)
        plt.plot(*lattice.singlets_to_plot(), "^", markersize=site_size)
        for num, coord in enumerate(zip(*lattice.oxygens_to_plot())):
            plt.annotate(str(num), coord)
        plt.axis("equal")
        plt.show()
        plt.plot(*lattice.points_to_plot(), "o", markersize=lattice_size)
        plt.plot(
            *lattice.midpoints_to_plot(),
            "x",
            markersize=site_size,
            markeredgewidth=cross_thickness,
        )
        plt.plot(*lattice.tripoints_to_plot(), "s", markersize=site_size)
        plt.plot(*lattice.singlets_to_plot(), "^", markersize=site_size)
        for num, coord in enumerate(zip(*lattice.points_to_plot())):
            plt.annotate(str(num), coord)
        plt.axis("equal")
        plt.show()
    elif plot:
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

    if similarity_filter:
        difference_distance = similarity_filter
    else:
        difference_distance = 0.0
    bit_lattice = lattice.get_intermediary(
        max_singlets=max_singlets,
        difference_distance=difference_distance,
        use_filter=similarity_filter is not None,
    )
    if rings_filter:
        noloops = lattice.no_rings()
        bit_lattice = bit_lattice.filtered(noloops)
    # print(bit_lattice)

    if use_parallel:
        solutions = bit_lattice.solve_parallel(True, silent=silent)
    else:
        solutions = bit_lattice.solve(True, silent=silent)

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
