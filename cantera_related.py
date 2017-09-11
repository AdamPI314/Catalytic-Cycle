"""
cantera related routines
"""
import os
import sys
import json
from collections import defaultdict
import numpy as np
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


def eval_2nd_order_rate_const(f_n, temp=None, pressure=1.0, buffer=None, rxn_idx=0):
    """
    f_n is path of "chem.xml" mechanism file
    temp is temperature
    buffer a dict store buffer gas and thier pressure
    rxn_idx is reaction index
    evaluate 2nd order rate costant
    """
    if buffer is None:
        buffer = dict()
    if temp is None:
        temp = [1000]

    spe = ct.Species.listFromFile(f_n)
    rxn = ct.Reaction.listFromFile(f_n)
    # print(S)
    rxn_str = []
    for r_t in rxn_idx:
        rxn_str.append(str(rxn[r_t]))
        # print(rxn[r_t])
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                      species=spe, reactions=rxn)

    # Units are a combination of kmol, m^3 and s, convert to cm^3/molecule/s
    # 10^6/(10^3 * 6.022 * 10^23)
    # kmol_m_s_to_molecule_cm_s = 10**6/(10**3 * 6.022 * 10**23)
    kmol_m_s_to_molecule_cm_s = 1e6 / (1e3 * 6.022e23)

    result = []
    for temp_t in temp:
        gas.TPY = temp_t, pressure * ct.one_atm, buffer
        # print(gas.forward_rates_of_progress[rxn_idx])
        result.append(
            gas.forward_rate_constants[rxn_idx] * kmol_m_s_to_molecule_cm_s)
    return rxn_str, result


def beta_1000_rate_constant_w2f(file_dir, beta=None, pressure=1.0, buffer=None, rxn_idx=0):
    """
    beta = 1000/T
    T = 1000/beta
    """
    if beta is None:
        beta = [1.0]
    np.savetxt(os.path.join(file_dir, "output", "beta.csv"),
               beta, fmt="%.15f", delimiter=",")
    temp = [1000.0 / b for b in beta]
    r_n, r_c = eval_2nd_order_rate_const(os.path.join(
        file_dir, "input", "chem.xml"), temp=temp, pressure=pressure, buffer=buffer, rxn_idx=rxn_idx)
    np.savetxt(os.path.join(file_dir, "output", "rxn_name.csv"),
               r_n, fmt="%s", delimiter=",")
    np.savetxt(os.path.join(file_dir, "output", "rate_constant.csv"),
               r_c, fmt="%1.15e", delimiter=",")


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # S_C = get_species_from_file(FILE_DIR)
    # write_spe_composition_2_file(FILE_DIR, S_C, tag="")
    BUFFER = {"npropyl": 1.0, "O2": 1.0, "HE": 1.0}
    RXN_IDX = [548, 549, 550, 551, 552, 553]
    # TEMP = [666.66666666]
    # BETA = [0.5, 1.0, 1.5, 2.0, 2.5]
    BETA = np.linspace(0.5, 2.5, num=25, endpoint=True)

    # eval_2nd_order_rate_const(os.path.join(
    #     FILE_DIR, "input", "chem.xml"), pressure=3.0, temp=TEMP, buffer=BUFFER, rxn_idx=RXN_IDX)
    beta_1000_rate_constant_w2f(
        FILE_DIR, beta=BETA, pressure=3.0, buffer=BUFFER, rxn_idx=RXN_IDX)
    print(FILE_DIR)
