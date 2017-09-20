"""
global settings
"""
import sys
import os
from collections import OrderedDict


def get_fast_rxn_trapped_spe(file_dir):
    """
    get_fast_rxn_trapped_spe
    """
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import fast_transitions
        return fast_transitions.get_fast_rxn_trapped_spe()
    except Exception:
        return OrderedDict(), OrderedDict()


def get_setting():
    """
    return global settings
    """
    setting = {
        # end time
        "end_t": 0.800000099855441071,
        "tau": 0.5,
        # atom followed
        "atom_f": "C",
        "init_s": 62,
        "end_s": "",
        # top n path
        "top_n_p": 1000,
        # top n species
        "top_n_s": 10,
        # number of trajectory used to generate pathway list running mc simulation
        "mc_n_traj": 10000000,
        # path integral number of trajectory
        "pi_n_traj": 10000,
        # number of time points when prepare path integral time points
        "pi_n_time": 1,
        # tag, M or fraction
        "tag": "M"
    }
    return setting
