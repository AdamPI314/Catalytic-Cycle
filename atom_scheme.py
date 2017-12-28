"""
followed atom scheme, either natural atoms, or hypothesized atoms
"""

import os
import sys
import time
from shutil import copy2
import read_write_configuration as rwc


def spe_composition_2_atom_scheme(file_dir):
    """
    convert species grouped atom scheme, which refers to file named
    "spe_composition.json" generated from cantera to a new file named
    "atom_scheme.json"
    """
    spe_comp = rwc.read_configuration(os.path.join(
        file_dir, "input", "spe_composition.json"))

    atom_scheme = {}
    for _, s_1 in enumerate(spe_comp):
        # print(s_1)
        for atom_1 in spe_comp[s_1]:
            if atom_1 not in atom_scheme:
                atom_scheme[atom_1] = {str(s_1): spe_comp[s_1][atom_1]}
            else:
                atom_scheme[atom_1].update({str(s_1): spe_comp[s_1][atom_1]})

    fn0 = os.path.join(file_dir, "input", "atom_scheme_backup.json")
    fn1 = os.path.join(file_dir, "input", "atom_scheme.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


def spe_information_2_atom_scheme(file_dir):
    """
    convert species information
    "species_information.json" to a new file named
    "atom_scheme.json"
    """
    spe_comp = rwc.read_configuration(os.path.join(
        file_dir, "input", "species_information.json"))

    atom_scheme = {}
    for _, spe_idx in enumerate(spe_comp):
        s_1 = spe_comp[spe_idx]["name"]
        for atom_1 in spe_comp[spe_idx]["spe_composition"]:
            atom_number = spe_comp[spe_idx]["spe_composition"][atom_1]
            if atom_number == "0":
                continue
            if atom_1 not in atom_scheme:
                atom_scheme[atom_1] = {str(s_1): float(atom_number)}
            else:
                atom_scheme[atom_1].update({str(s_1): float(atom_number)})

    fn0 = os.path.join(file_dir, "input", "atom_scheme_backup.json")
    fn1 = os.path.join(file_dir, "input", "atom_scheme.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(atom_scheme, fn1)


if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    # spe_composition_2_atom_scheme(FILE_DIR)
    spe_information_2_atom_scheme(FILE_DIR)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
