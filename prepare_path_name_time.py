"""
prepare pathway name and pathway time for pathway-probability evaluation
"""
import os
import sys
import numpy as np
import pandas as pd
from copy import deepcopy
import parse_spe_reaction_info as psri
import parse_pattern as pp


def prepare_pathway_name(
        data_dir, top_n=5, flag="", delimiter=",", end_s_idx=None, species_path=False,
        path_reg=None, spe_idx=None):
    """
    prepare pathway_name_candidate.csv
    """
    # read from pathway_stat.csv
    prefix = ""
    if species_path is True:
        prefix = "species_"
    f_n_ps = os.path.join(data_dir, "output", prefix + "pathway_stat.csv")

    if flag == "":
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate.csv")
    else:
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate_" + str(flag) + ".csv")

    try:
        os.remove(f_n_pn)
    except OSError:
        pass

    path_list = []
    d_f = pd.read_csv(f_n_ps, names=['pathway', 'frequency'])
    if path_reg is not None:
        mask1 = d_f['pathway'].str.contains(path_reg)
    else:
        mask1 = d_f['pathway'].str.len() > 0

    if spe_idx is None:
        mask2 = d_f['pathway'].str.len() > 0
    else:
        net_product = psri.parse_reaction_net_product(data_dir)
        s_p_r_c = psri.parse_species_pair_reaction(data_dir)

        def func_filter(x): return pp.parse_species_along_path_using_reaction(
            pathname=x, net_r_p=net_product, spe_idx=spe_idx, s_p_r_c=s_p_r_c) >= 1
        mask2 = d_f.apply(func_filter, axis=1)

    # read
    if end_s_idx is None or end_s_idx == []:
        mask3 = d_f['pathway'].str.len() > 0
        path_list.extend(d_f[mask1 & mask2 & mask3]['pathway'][0:top_n])
    else:
        for s_i in end_s_idx:
            mask3 = d_f['pathway'].str.endswith("S" + str(s_i))
            path_list.extend(d_f[mask1 & mask2 & mask3]['pathway'][0:top_n])

    # save
    np.savetxt(f_n_pn, path_list, fmt="%s")
    return len(path_list)


def prepare_pathway_name_for_passage_time(data_dir, flag="", init_s_idx=None):
    """
    prepare pathway_name_candidate.csv
    """
    # read from pathway_stat.csv
    prefix = "species_"

    if flag == "":
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate.csv")
    else:
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate_" + str(flag) + ".csv")

    try:
        os.remove(f_n_pn)
    except OSError:
        pass

    # read
    if init_s_idx is None:
        init_s_idx_tmp = [62]
    else:
        init_s_idx_tmp = deepcopy(init_s_idx)

    path_list = []
    for s_i in init_s_idx_tmp:
        path_list.extend(["S" + str(s_i) + "R100S100"])

    # save
    np.savetxt(f_n_pn, path_list, fmt="%s")
    return len(path_list)


def prepare_pathway_time(data_dir, top_n=5, num=1, flag="", begin_t=0.0, end_t=1.0, species_path=False):
    """
    prepare pathway_time.csv
    num represents number of points
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    if flag == "":
        f_n_pt = os.path.join(data_dir, "output",
                              prefix + "pathway_time_candidate.csv")
    else:
        f_n_pt = os.path.join(data_dir, "output",
                              prefix + "pathway_time_candidate_" + str(flag) + ".csv")

    try:
        os.remove(f_n_pt)
    except OSError:
        pass

    # time matrix
    t_mat = np.empty((top_n, num + 1, ))
    for idx, _ in enumerate(t_mat):
        t_mat[idx] = np.linspace(begin_t, end_t, num + 1)

    np.savetxt(f_n_pt, t_mat[:, 1::], delimiter=',', fmt='%.15f')


if __name__ == '__main__':
    #     print("hello")
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
#     print(DATA_DIR)

    # prepare_pathway_name(DATA_DIR, top_n=5, flag="",
    #                      delimiter=",", end_s_idx=[62, 59])
    prepare_pathway_name(DATA_DIR, top_n=10, flag="", delimiter=',',
                         end_s_idx=None, species_path=False, path_reg='^S62R[736|738]')
