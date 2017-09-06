"""
cantera related routines
"""
import os
import sys
import json
from collections import defaultdict
try:
    sys.path.append('/opt/cantera/v2.4.0/lib/python3.5/site-packages')
    import cantera as ct
except ImportError:
    print("CANTERA NO AVAILABLE")


def get_species_from_file(file_dir):
    """
    get species composition
    """
    all_species = ct.Species.listFromFile(
        os.path.join(file_dir, "input", "chem.xml"))
    spe_composition = defaultdict(dict)

    for spe in all_species:
        spe_composition[spe.name] = spe.composition
    return spe_composition


def write_spe_composition_2_file(file_dir, spe_composition, tag=None):
    """
    write to json file
    """
    if tag is not None and tag != "":
        f_n = "spe_composition" + str(tag) + ".json"
    else:
        f_n = "spe_composition.json"
    with open(os.path.join(file_dir, "input", f_n), 'w') as f_h:
        json.dump(spe_composition, f_h, indent=4, sort_keys=False)


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    S_C = get_species_from_file(FILE_DIR)
    write_spe_composition_2_file(FILE_DIR, S_C, tag="")
    print(FILE_DIR)
