import numpy as np
import pandas as pd
import sys
import re
import os

import parse_spe_reaction_info as psri

if __name__ == '__main__':
    print("test")
    file_dir = os.path.abspath(os.path.realpath(
        os.path.join(sys.argv[0], os.pardir, os.pardir, os.pardir)))
    print(file_dir)

    filename_p = os.path.join(file_dir, "output", "pathway_stat.csv")
    # load data
    path_data = pd.read_csv(filename_p, names=['path', 'prob'])
    print(path_data.head())
    # renormalize path probability
    total_prob = sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    path_data['prob'] /= total_prob
    print(path_data.head())
    # parse species and reaction info

    filename1 = os.path.join(file_dir, "output", "species_labelling.csv")
    filename2 = os.path.join(file_dir, "output", "reaction_labelling.csv")
    spe_ind_name_dict, spe_name_ind_dict = psri.parse_spe_info(filename1)
    old_new_index_dict, new_old_index_dict, new_ind_reaction_dict = psri.parse_reaction_and_its_index(
        filename2)

    print(spe_ind_name_dict)
    print(new_ind_reaction_dict)
    # convert path name to real species and reaction instead of index
    # use pandas map function
    path_data['path'] = path_data['path'].apply(
        lambda x: psri.PATH_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, x))
