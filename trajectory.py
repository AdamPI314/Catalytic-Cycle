"""
get trajectory related information
"""

import os
import sys
import numpy as np
from collections import defaultdict, OrderedDict
import parse_spe_reaction_info as psri


def get_species_concentration(file_dir, top_n=10, tau=1.0):
    """
    get species concentration at a tau, where tau is the ratio of the time_wanted/end_time
    """
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_fraction.csv"), delimiter=",")
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

    for idx, val in enumerate(c_idx_map):
        if idx < top_n:
            spe_idx = next(iter(c_idx_map[val]))
            print(idx, val, spe_idx, spe_idx_name_dict[spe_idx])


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    get_species_concentration(FILE_DIR, top_n=20, tau=0.9)
