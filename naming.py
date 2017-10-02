"""
namespace, naming stuff
"""

import os
import sys
from shutil import copyfile
import read_write_configuration as rwc


def get_suffix(file_dir, init_spe=None, atom_followed=None, tau=None, pathwayEndWith=None):
    """
    get suffix
    """
    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))
    suffix = ""
    if init_spe is None:
        suffix += "_S" + str(setting['pathway']['init_spe'])
    else:
        suffix += "_S" + str(init_spe)
    if atom_followed is None:
        suffix += "_" + str(setting['pathway']['atom_followed'])
    else:
        suffix += "_" + str(atom_followed)
    if tau is None:
        suffix += "_" + str(setting['pathway']['tau'])
    else:
        suffix += "_" + str(tau)
    if pathwayEndWith is None:
        suffix += "_" + str(setting['pathway']['pathwayEndWith'])
    else:
        suffix += "_" + str(pathwayEndWith)

    return suffix


def copy_sohr_files(file_dir):
    """
    make a copy of SOHR files
    1. output/pathway_stat.csv
    2. input/pathway_name.csv
    3. input/pathway_time.csv
    4. output/pathway_name.csv
    5. output/pathway_prob.csv
    6. output/pathname_prob.csv
    7. output/chattering_group_info.json
    """
    suffix = get_suffix(file_dir)

    f_n_1 = os.path.join(file_dir, "output", "pathway_stat.csv")
    f_n_2 = os.path.join(file_dir, "output", "pathway_stat" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_2 = os.path.join(file_dir, "output", "pathway_name" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", "pathway_prob.csv")
    f_n_2 = os.path.join(file_dir, "output", "pathway_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", "pathname_prob.csv")
    f_n_2 = os.path.join(file_dir, "output", "pathname_prob" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "output", "chattering_group_info.json")
    f_n_2 = os.path.join(file_dir, "output", "chattering_group_info" + suffix + ".json")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "input", "pathway_name.csv")
    f_n_2 = os.path.join(file_dir, "input", "pathway_name" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    f_n_1 = os.path.join(file_dir, "input", "pathway_time.csv")
    f_n_2 = os.path.join(file_dir, "input", "pathway_time" + suffix + ".csv")
    if os.path.isfile(f_n_1):
        print(f_n_1, "found")
        copyfile(f_n_1, f_n_2)

    return


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    copy_sohr_files(FILE_DIR)
