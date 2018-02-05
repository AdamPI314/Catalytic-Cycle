"""
to deal with chattering
"""

import sys
import os
import time
from shutil import copy2
import numpy as np
import parse_spe_reaction_info as psri
import read_write_configuration as rwc
import interpolation


def initiate_fast_reaction(file_dir):
    """
    intiate file named "reaction_info_base.json" based on file named "reaction_labelling.csv"
    """
    new_old_index_dict, new_ind_reaction_dict = psri.parse_reaction_and_its_index(
        file_dir)

    rxn_pair_dict = dict()

    # un-paired species
    unpaired = dict()
    for _, val1 in enumerate(new_old_index_dict):
        # print(idx, val1, new_old_index_dict[val1])
        this_value = int(new_old_index_dict[val1])
        neg_value = -1 * this_value
        if (neg_value in unpaired):
            # print(int(val1), unpaired[neg_value])
            rxn_pair_dict.update({int(val1): unpaired[neg_value]})
            rxn_pair_dict.update({unpaired[neg_value]: int(val1)})
            unpaired.pop(neg_value)
        else:
            unpaired.update({this_value: int(val1)})

    rxn_info = {}
    for _, val1 in enumerate(new_ind_reaction_dict):
        entry = {
            str(val1): {
                "formula": new_ind_reaction_dict[val1],
                # default value, 10^-100, assume the reaction is super slow
                "time_scale": -100,
                "reverse_reaction": "None"
            }
        }
        if int(val1) in rxn_pair_dict:
            entry[str(val1)]["reverse_reaction"] = \
                str(rxn_pair_dict[int(val1)])
        rxn_info.update(entry)

    fn0 = os.path.join(file_dir, "input", "reaction_info_base_backup.json")
    fn1 = os.path.join(file_dir, "input", "reaction_info_base.json")

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(rxn_info, fn1)


def update_fast_reaction(file_dir, tau=0.7, end_t=1.0, tag="M"):
    """
    update fast reaction based on reference trajectory
    """
    fn0 = os.path.join(file_dir, "input", "reaction_info_base_backup.json")
    fn1 = os.path.join(file_dir, "input", "reaction_info_base.json")

    rxn_info = rwc.read_configuration(fn1)

    time_v = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    rxn_rates = np.loadtxt(os.path.join(file_dir, "output",
                                        "reaction_rate_dlsode_" + str(tag) + ".csv"), delimiter=",")

    actual_time = float(tau) * float(end_t)
    for _, val in enumerate(rxn_info):
        actual_rate = interpolation.interp1d(
            time_v, rxn_rates[:, int(val)], actual_time)
        if actual_rate != 0:
            time_scale = np.log10(actual_rate)
            if time_scale >= -100 and time_scale <= 100:
                # print(time_scale)
                rxn_info[val]["time_scale"] = time_scale

    if os.path.isfile(fn1):
        copy2(fn1, fn0)

    rwc.write_configuration(rxn_info, fn1)


def fast_reaction_w2f(file_dir, threshold=-7):
    """
    prepare "fast_reaction.json" file, this file will be
    1) manually changed later
    2) used to generate chattering information later
    """
    fast_transition = {}
    counter = 0

    fn_rib = os.path.join(file_dir, "input", "reaction_info_base.json")
    rxn_info = rwc.read_configuration(fn_rib)

    # care fast reaction pair only, both forward and backward reactions are fast
    unpaired_fast_rxn = set()
    for _, val in enumerate(rxn_info):
        if rxn_info[val]["reverse_reaction"] != "None" \
                and float(rxn_info[val]["time_scale"]) >= threshold:
            # print(rxn_info[val]["formula"])
            this_rxn = str(val)
            paired_rxn = rxn_info[val]["reverse_reaction"]
            if paired_rxn not in unpaired_fast_rxn:
                unpaired_fast_rxn.add(this_rxn)
            else:
                entry = {
                    str(counter): {
                        "formula1": rxn_info[paired_rxn]["formula"],
                        "formula2": rxn_info[this_rxn]["formula"],
                        "reaction1": int(paired_rxn),
                        "reaction2": int(this_rxn)
                    }
                }
                fast_transition.update(entry)
                counter += 1

    fn_frb0 = os.path.join(file_dir, "input", "fast_reaction_base_backup.json")
    fn_frb1 = os.path.join(file_dir, "input", "fast_reaction_base.json")
    if os.path.isfile(fn_frb1):
        copy2(fn_frb1, fn_frb0)

    rwc.write_configuration(fast_transition, fn_frb1)


def generate_fast_rxn_chattering_spe(file_dir):
    """
    generate fast reaction and chattering species based on four files
    0) species_information.json
    1) reaction_information.json
    2) atom_scheme.json
    3) fast_reaction_base.json

    save file named "fast_transition.json", this file will be used to update file
    "local_settings.py" manually
    """
    f_n_si = os.path.join(file_dir, "input", "species_information.json")
    f_n_ri = os.path.join(file_dir, "input", "reaction_information.json")
    f_n_as = os.path.join(file_dir, "input", "atom_scheme.json")
    f_n_frb = os.path.join(file_dir, "input", "fast_reaction_base.json")

    spe_info = rwc.read_configuration(f_n_si)
    rxn_info = rwc.read_configuration(f_n_ri)
    atom_scheme = rwc.read_configuration(f_n_as)
    fast_rxn_base = rwc.read_configuration(f_n_frb)

    fast_transition = []

    for _, val1 in enumerate(fast_rxn_base):
        entry = {}

        rxn_1_idx = fast_rxn_base[val1]["reaction1"]
        rxn_2_idx = fast_rxn_base[val1]["reaction2"]

        reactant1 = rxn_info[str(rxn_1_idx)]["net_reactant"]
        reactant2 = rxn_info[str(rxn_2_idx)]["net_reactant"]

        entry.update({"formula1": fast_rxn_base[val1]["formula1"]})
        entry.update({"formula2": fast_rxn_base[val1]["formula2"]})
        entry.update({"rxn": [int(rxn_1_idx), int(rxn_2_idx)]})
        entry.update({"spe": {}})

        s_1_idx = "None"
        s_2_idx = "None"
        for atom_followed in atom_scheme:
            for _, val2 in enumerate(reactant1):
                spe_idx = reactant1[val2]["species_index"]
                spe_name = spe_info[spe_idx]["name"]
                if spe_name in atom_scheme[atom_followed]:
                    s_1_idx = str(spe_idx)
                    # only one species allowed
                    break
            for _, val2 in enumerate(reactant2):
                spe_idx = reactant2[val2]["species_index"]
                spe_name = spe_info[spe_idx]["name"]
                if spe_name in atom_scheme[atom_followed]:
                    s_2_idx = str(spe_idx)
                    # only one species allowed
                    break
            if s_1_idx != "None" and s_2_idx != "None":
                entry["spe"].update(
                    {atom_followed: [int(s_1_idx), int(s_2_idx)]})

        fast_transition.append(entry)

    fn_ft0 = os.path.join(file_dir, "input", "fast_transition_backup.json")
    fn_ft1 = os.path.join(file_dir, "input", "fast_transition.json")
    if os.path.isfile(fn_ft1):
        copy2(fn_ft1, fn_ft0)

    rwc.write_configuration(fast_transition, fn_ft1)


if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(FILE_DIR)

    # initiate_fast_reaction(FILE_DIR)
    # update_fast_reaction(FILE_DIR, tau=0.7, end_t=0.25)
    # fast_reaction_w2f(FILE_DIR, threshold=-8)
    generate_fast_rxn_chattering_spe(FILE_DIR)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
