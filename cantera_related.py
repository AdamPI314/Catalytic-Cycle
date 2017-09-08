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
    print("CANTERA NOT AVAILABLE")


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


def evaluate_rate_constant(f_n, temp=1000, buffer=None, rxn_idx=0):
    """
    f_n is path of "chem.xml" mechanism file
    temp is temperature
    buffer a dict store buffer gas and thier pressure
    rxn_idx is reaction index
    """
    if buffer is None:
        buffer = dict()

    S = ct.Species.listFromFile(f_n)
    R = ct.Reaction.listFromFile(f_n)
    print(S)
    print(R)
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                      species=S, reactions=R)

    gas.TPY = 650, 3 * ct.one_atm, 'npropyl:1.0, O2:1.0, HE:1.0'

    print(gas.reactant_stoich_coeffs())
    print(gas.product_stoich_coeffs())

    print(gas.reverse_rate_constants)
    print(gas.forward_rates_of_progress)

    return


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # S_C = get_species_from_file(FILE_DIR)
    # write_spe_composition_2_file(FILE_DIR, S_C, tag="")
    BUFFER = {"npropyl": 1.0, "O2": 1.0, "HE": 1.0}
    evaluate_rate_constant(os.path.join(
        FILE_DIR, "input", "chem.xml"), temp=650, buffer=BUFFER, rxn_idx=0)
    print(FILE_DIR)
