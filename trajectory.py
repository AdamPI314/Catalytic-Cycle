"""
get trajectory related information
"""

import os
import sys
from collections import defaultdict, OrderedDict
import numpy as np
import parse_spe_reaction_info as psri


def get_species_with_top_n_concentration(file_dir, exclude, top_n=10, tau=1.0, tag="fraction", atoms=None):
    """
    get species concentration at a tau, where tau is the ratio of the time_wanted/end_time
    """
    if atoms is None:
        atoms = ["C"]
    if exclude is None:
        exclude = []
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")
    time_idx = int(tau * len(conc))
    if time_idx >= len(conc):
        time_idx = (len(conc) - 1)

    data = conc[time_idx, :]
    c_idx_map = defaultdict(set)
    for idx, val in enumerate(data):
        c_idx_map[val].add(str(idx))
    c_idx_map = OrderedDict(sorted(c_idx_map.items(), reverse=True))

    spe_idx_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    spe_composition = psri.read_spe_composition(
        os.path.join(file_dir, "input", "spe_composition.json"))

    spe_idx_list = []
    counter = 0

    for _, val in enumerate(c_idx_map):
        if counter < top_n:
            spe_idx = next(iter(c_idx_map[val]))
            indicator = False
            for _, atom in enumerate(atoms):
                if atom in spe_composition[spe_idx_name_dict[spe_idx]]:
                    indicator = True
                    break
            if spe_idx_name_dict[spe_idx] not in exclude and indicator:
                # print(val, spe_idx, spe_idx_name_dict[spe_idx])
                spe_idx_list.append(int(spe_idx))
                counter += 1

    # species doesn't contain atom we are interested in
    exclude_spe_name_list = []
    for idx, s_n_t in enumerate(spe_composition):
        indicator = False
        for _, atom in enumerate(atoms):
            if atom in spe_composition[s_n_t]:
                indicator = True
        if indicator is False:
            exclude_spe_name_list.append(s_n_t)

    return spe_idx_list, exclude_spe_name_list


def get_normalized_concentration(file_dir, tag="fraction", exclude_names=None, renormalization=True):
    """
    return normalized concentration
    """
    if exclude_names is None:
        exclude_names = []
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")
    if renormalization is False:
        return conc

    _, s_n_idx = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    # renormalization
    exclude_idx_list = [int(s_n_idx[x]) for x in exclude_names]
    # set the concentration of these species to be zero
    for _, idx in enumerate(exclude_idx_list):
        conc[:, idx] = 0.0
    # normalize
    for idx, _ in enumerate(conc):
        conc[idx, :] /= np.sum(conc[idx, :])

    return conc


def get_normalized_concentration_at_time(file_dir, tag="fraction", tau=1.0, exclude_names=None, renormalization=True):
    """
    return normalized species concentration at time
    """
    if exclude_names is None:
        exclude_names = []

    f_n = os.path.join(file_dir, "output", "spe_concentration_" +
                       str(tau) + "_dlsode_" + str(tag) + ".csv")
    conc = np.loadtxt(f_n, delimiter=",")

    if renormalization is False:
        return conc

    _, s_n_idx = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    # renormalization
    exclude_idx_list = [int(s_n_idx[x]) for x in exclude_names]
    # set the concentration of these species to be zero
    for _, idx in enumerate(exclude_idx_list):
        conc[idx] = 0.0
    conc /= np.sum(conc)
    return conc


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    # get_species_with_top_n_concentration(
    #     FILE_DIR, exclude=["N2", "AR"], top_n=10, tau=0.9, tag="M", atoms=["C"])
    get_normalized_concentration_at_time(
        FILE_DIR, tag="M", tau=0.9, exclude_names=None, renormalization=True)
