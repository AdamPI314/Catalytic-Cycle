"""
to deal with chattering
"""

import sys
import os
import time
import parse_spe_reaction_info as psri


def initiate_fast_reaction(file_dir):
    """
    intiate file named "reaction_info.json" based on file named "reaction_labelling.csv"
    """
    new_old_index_dict, new_ind_reaction_dict = psri.parse_reaction_and_its_index(
        os.path.join(file_dir, "input", "reaction_labelling.csv"))
    print(new_old_index_dict, new_ind_reaction_dict)


if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    initiate_fast_reaction(FILE_DIR)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
