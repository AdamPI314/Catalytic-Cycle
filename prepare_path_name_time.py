"""
prepare pathway name and pathway time for pathway-probability evaluation
"""
import os
import sys
from copy import deepcopy
import numpy as np
import pandas as pd
import parse_spe_reaction_info as psri
import parse_pattern as pp


def prepare_pathway_name(
        data_dir, top_n=5, flag="", end_s_idx=None, species_path=False,
        path_reg=None, spe_idx=None, spe_production_oriented=False):
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
        net_reactant = psri.parse_reaction_net_reactant(data_dir)
        net_product = psri.parse_reaction_net_product(data_dir)
        s_p_r_c = psri.parse_species_pair_reaction(data_dir)

        mask2 = d_f.apply(lambda x: pp.parse_net_species_along_path_using_reaction(
            pathname=x['pathway'], net_r=net_reactant, net_p=net_product, spe_idx=spe_idx, s_p_r_c=s_p_r_c) >= 1, axis=1)

    # read
    if end_s_idx is None or end_s_idx == []:
        mask3 = d_f['pathway'].str.len() > 0
        path_list.extend(d_f[mask1 & mask2 & mask3]['pathway'])

        if spe_production_oriented is False:
            # save
            np.savetxt(f_n_pn, path_list[0:top_n], fmt="%s")
            return len(path_list[0:top_n])
        elif spe_idx is not None:
            path_list2 = []
            path_set = set()

            for _, val1 in enumerate(path_list):
                p_list, r_list = pp.get_spe_production_sub_path(
                    val1, net_product, spe_idx, s_p_r_c)
                for idx2, val2 in enumerate(r_list):
                    if val2 not in path_set:
                        path_set.add(val2)
                        path_list2.append(p_list[idx2])

            np.savetxt(f_n_pn, path_list2[0:top_n], fmt="%s")
            return len(path_list2[0:top_n])
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


def prepare_pathway_time(
        data_dir, top_n=5, num=1, flag="", begin_t=0.0, end_t=1.0,
        species_path=False, fixed_t0_or_tf=None):
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

    if fixed_t0_or_tf is None or fixed_t0_or_tf == "t0":
        np.savetxt(f_n_pt, t_mat[:, 1::], delimiter=',', fmt='%.15f')
    else:
        np.savetxt(f_n_pt, t_mat[:, :-1], delimiter=',', fmt='%.15f')


if __name__ == '__main__':
    #     print("hello")
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
#     print(DATA_DIR)

    # prepare_pathway_name(DATA_DIR, top_n=5, flag="",
    #                      end_s_idx=[62, 59])
    # prepare_pathway_name(DATA_DIR, top_n=10, flag="",
    #                      end_s_idx=None, species_path=False, path_reg='^S62R[736|738]')

    prepare_pathway_name(
        DATA_DIR, top_n=5, flag="", end_s_idx=None, species_path=False,
        path_reg=None, spe_idx=10, spe_production_oriented=True)
