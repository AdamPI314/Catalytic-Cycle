"""
namespace, naming stuff
"""

import os
import sys
from shutil import copyfile
import global_settings


def valid_name(name1):
    name2 = ''
    for c in name1:
        if (c >= 'A' and c <= 'Z') or (c >= 'a' and c <= 'z') or (c >= '0' and c <= '9'):
            name2 += c
        elif len(name2) > 0 and name2[-1] != '_':
            name2 += '_'
    # get rid of trailing '_'
    if len(name2) > 0 and name2[-1] == '_':
        return name2[0:-1]
    return name2


def get_suffix(data_dir, init_spe=None, atom_followed=None, end_t=None):
    """
    get suffix
    """
    g_s = global_settings.get_setting(data_dir)
    suffix = ""
    if init_spe is None:
        suffix += "_S" + str(g_s['init_s'])
    else:
        suffix += "_S" + str(init_spe)
    if atom_followed is None:
        suffix += "_" + str(g_s['atom_f'])
    else:
        suffix += "_" + str(atom_followed)
    if end_t is None:
        suffix += "_" + str(g_s['end_t'])
    else:
        suffix += "_" + str(end_t)
    # add mc_n_traj
    suffix += "_" + str(int(g_s['mc_n_traj']))
    # add pi_n_traj
    suffix += "_" + str(int(g_s['pi_n_traj']))
    # add top_n_p
    suffix += "_" + str(int(g_s['top_n_p']))
    # path reg
    if g_s['path_reg'] is not None:
        suffix += "_include_" + str(valid_name(g_s['path_reg']))
    # no path reg
    if g_s['no_path_reg'] is not None:
        suffix += "_exclude_" + str(valid_name(g_s['no_path_reg']))

    return suffix


def copy_sohr_files(data_dir, species_path=False):
    """
    make a copy of SOHR files
    1. output/pathway_stat.csv
    2. output/pathway_name_candidate.csv
    3. output/pathway_time_candidate.csv
    4. output/pathway_name_candidate_real_path.csv
    5. output/pathway_prob.csv
    6. output/species_pathway_AT.csv
    7. output/pathname_prob.csv
    8. output/chattering_group_info.json
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = get_suffix(data_dir)

    f_n_1 = os.path.join(data_dir, "output", prefix + "pathway_stat.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "pathway_stat" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix + "pathway_prob.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "pathway_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix + "species_pathway_AT.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "species_pathway_AT" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "species_pathway_AT_no_IT.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "species_pathway_AT_no_IT" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "species_pathway_AT_with_SP.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "species_pathway_AT_with_SP" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix + "pathway_SP.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "pathway_SP" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix + "pathname_prob.csv")
    f_n_2 = os.path.join(data_dir, "output", prefix +
                         "pathname_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", "chattering_group_info.json")
    f_n_2 = os.path.join(data_dir, "output",
                         "chattering_group_info" + suffix + ".json")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "pathway_name_candidate.csv")
    f_n_2 = os.path.join(data_dir, "output",
                         prefix + "pathway_name_candidate" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "pathway_time_candidate.csv")
    f_n_2 = os.path.join(data_dir, "output",
                         prefix + "pathway_time_candidate" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "pathway_name_candidate_real_path.csv")
    f_n_2 = os.path.join(data_dir, "output",
                         prefix + "pathway_name_candidate_real_path" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(data_dir, "output", prefix +
                         "Merchant_f_2d.csv")
    f_n_2 = os.path.join(data_dir, "output",
                         prefix + "Merchant_f_2d" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)
    return


if __name__ == '__main__':
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)
    # copy_sohr_files(DATA_DIR)
    print(get_suffix(DATA_DIR))
