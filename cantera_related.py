"""
cantera related routines
"""
import os
import sys
import json
from collections import defaultdict
try:
    import cantera as ct
except ImportError:
    print("CANTERA NOT AVAILABLE")


def get_species_from_file(file_dir):
    """
    get species composition
    """
    all_species = ct.Species.listFromFile(
        os.path.join(file_dir, "input", "chem.cti"))
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


def evaluate_rate_constant(f_n, temp=1000, pressure=1.0, buffer=None, rxn_idx=0):
    """
    f_n is path of "chem.xml" mechanism file
    temp is temperature
    buffer a dict store buffer gas and thier pressure
    rxn_idx is reaction index
    """
    if buffer is None:
        buffer = dict()

    spe = ct.Species.listFromFile(f_n)
    rxn = ct.Reaction.listFromFile(f_n)
    # print(S)
    for r in rxn_idx:
        print(rxn[r])
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                      species=spe, reactions=rxn)

    # gas.TPY = temp, pressure * ct.one_atm, 'npropyl:1.0, O2:1.0, HE:1.0'
    gas.TPY = temp, pressure * ct.one_atm, buffer

    # print(gas.reactant_stoich_coeffs())
    # print(gas.product_stoich_coeffs())

    # Units are a combination of kmol, m^3 and s, convert to cm^3/molecule/s
    # 10^6/(10^3 * 6.022 * 10^23)
    # kmol_m_s_to_molecule_cm_s = 10**6/(10**3 * 6.022 * 10**23) 
    kmol_m_s_to_molecule_cm_s = 1e6/(1e3 * 6.022e23) 
    print(gas.forward_rate_constants[rxn_idx] * kmol_m_s_to_molecule_cm_s)
    # print(gas.reverse_rate_constants)
    # print(gas.forward_rates_of_progress[rxn_idx])

    return


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # S_C = get_species_from_file(FILE_DIR)
    # write_spe_composition_2_file(FILE_DIR, S_C, tag="")
    BUFFER = {"npropyl": 1.0, "O2": 1.0, "HE": 1.0}
    RXN_IDX = [548, 549, 550, 551, 552, 553]
    TEMP = 666.66666666
    evaluate_rate_constant(os.path.join(
        FILE_DIR, "input", "chem.xml"), pressure=3.0, temp=TEMP, buffer=BUFFER, rxn_idx=RXN_IDX)
    print(FILE_DIR)
