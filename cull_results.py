import os
import json
from itertools import chain, combinations
from math import sqrt
from alive_progress import alive_bar


# FOCUS = {"sizeable_rust_32232.json", "sizeable_rust_35469.json"}
FOCUS = set()


def nearly_in(new: list[float], uniques: list[list[float]], margin: float) -> bool:
    """
    Returns true if `new` list has matching list in `uniques`, where all elements 
    differ by no more than `margin`:

    For every `x`:
    ```
    new[x] - uniques[:][x] <= margin
    ```
    """
    for unique in uniques:
        if not len(new) == len(unique):
            continue
        for a, b in zip(new, unique):
            if (a - b) > margin:
                break
        else:
            return True
    return False


def sorting_key(input) -> float:
    return 100 * input["x"] + input["y"]


def cull(path: str):
    list_of_files = os.listdir(path)
    new_path = path + "_unique"
    os.mkdir(new_path)
    unique_spacings = {}
    unique_files = {}
    total_amount = 0
    
    with open(f"{path}/{list_of_files[0]}") as file:
        lattice = json.load(file)
        max_x = 0
        max_y = 0
        for number, lattice_point in enumerate(lattice["lattice_points"]):
            if lattice_point["x"] > max_x and lattice_point["y"] == 0:
                print(f"Choose point {number} for {max_x=}")
                max_x = lattice_point["x"]
            if lattice_point["y"] > max_y:
                print(f"Choose point {number} for {max_y=}")
                max_y = lattice_point["y"]
    # extremely important progress bar
    with alive_bar(
        len(list_of_files), spinner="dots_waves", title="Analyzing files..."
    ) as bar:
        for filename in list_of_files:
            with open(f"{path}/{filename}") as file:
                lattice = json.load(file)

                connection_type_count = (
                    len(lattice["tripoints"]),
                    len(lattice["midpoints"]),
                    len(lattice["singles"]),
                )
                oxygens = list(
                    chain(
                        lattice["tripoints"], lattice["midpoints"], lattice["singles"]
                    )
                )
                
                oxygen_spacing = []
                for o1, o2 in combinations(oxygens, 2):
                    dx = abs(o1["x"] - o2["x"])
                    if dx > max_x/2:
                        dx -= max_x
                    dy = abs(o1["y"] - o2["y"])
                    if dy > max_y/2:
                        dy -= max_y

                    spacings = sqrt((dx)**2 + (dy)**2)
                    oxygen_spacing.append(spacings)
                oxygen_spacing.sort()
                if filename in FOCUS:
                    print([round(number, 4) for number in oxygen_spacing])

                if not unique_spacings.get(connection_type_count):
                    unique_spacings[connection_type_count] = [oxygen_spacing]
                    unique_files[connection_type_count] = [filename]
                    total_amount += 1
                elif not nearly_in(
                    oxygen_spacing, unique_spacings[connection_type_count], 0.0001
                ):
                    unique_spacings[connection_type_count].append(oxygen_spacing)
                    unique_files[connection_type_count].append(filename)
                    total_amount += 1
            bar()

    # print(total_amount)
    # print(f"{len(unique_files[(4, 0, 4)])}")
    with alive_bar(
        total_amount, spinner="dots_waves", title="Copying uniques..."
    ) as bar:
        for conn_type, files in unique_files.items():
            for filename in files:
                if f"{conn_type}" not in os.listdir(new_path):
                    os.mkdir(f"{new_path}/{conn_type}")
                with open(f"{path}/{filename}") as file:
                    with open(f"{new_path}/{conn_type}/{filename}", "w") as new_file:
                        json.dump(json.load(file), new_file)
                bar()


if __name__ == "__main__":
    cull("exports/sizeable")
