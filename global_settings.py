"""
global settings
"""


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
        "top_n_p": 5000,
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
