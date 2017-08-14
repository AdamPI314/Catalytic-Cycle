"""
parse species patten, find catalytic cycle
"""

import os
import sys
import re
from itertools import combinations
import numpy as np
import pandas as pd
import parse_spe_reaction_info as psri


def parse_species(path):
    """
    parse sepcies name, return a dictionary of species and count
    """
    matched_tmp = re.findall(r"(S\d+)", path)
    d_map = dict()
    # not the first and last element
    for _, val in enumerate(matched_tmp[1:-1]):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_reaction(path):
    """
    parse reactions, return a dictionary of reactions and their count
    """
    matched_tmp = re.findall(r"(R\d+)", path)
    d_map = dict()
    # not the first and last element
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_species_cycle(path):
    """
    parse sepcies name, return a dictionary of species cycle and count
    """
    matched_tmp = re.findall(r"(S\d+)", path)
    d_map = dict()
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    cycle_map = dict()
    for key, val in d_map.items():
        # search pathway if occurence of spcies >= 2
        if val >= 2:
            idx_tmp = [(m.start(0), m.end(0)) for m in re.finditer(key, path)]
            for p_x in combinations(idx_tmp, 2):
                cycle = path[p_x[0][0]:p_x[1][1]]
                if cycle not in cycle_map:
                    cycle_map[cycle] = 1
                else:
                    cycle_map[cycle] += 1

    return cycle_map


def species_count(file_dir, top_n=50):
    """
    species occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_species(p_n)
        for key, value in map_tmp.items():
            if key not in spe_map:
                spe_map[key] = value * p_p
            else:
                spe_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(spe_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'species', 'frequency'])
    f_n_out1 = os.path.join(file_dir, "output", "species_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])

    # load spe and reaction info
    spe_ind_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

    # convert species reaction index to real species and reactions
    d_f['species'] = d_f['species'].apply(
        lambda x: psri.pathname_to_real_spe_reaction(
            spe_ind_name_dict, new_ind_reaction_dict, x)
        .strip())
    f_n_out2 = os.path.join(file_dir, "output", "species_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])


def reaction_count(file_dir, top_n=50):
    """
    reaction occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_reaction(p_n)
        for key, value in map_tmp.items():
            if key not in spe_map:
                spe_map[key] = value * p_p
            else:
                spe_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(spe_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'reaction', 'frequency'])
    f_n_out1 = os.path.join(file_dir, "output", "reaction_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    print(d_f['reaction'])
    f_n_out2 = os.path.join(file_dir, "output", "reaction_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


def species_cycle(file_dir, top_n=50):
    """
    species cycle in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_cycle_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_species_cycle(p_n)
        for key, value in map_tmp.items():
            if key not in spe_cycle_map:
                spe_cycle_map[key] = value * p_p
            else:
                spe_cycle_map[key] += value * p_p

    d_f = pd.DataFrame(list(sorted(spe_cycle_map.items(), key=lambda x: x[1], reverse=True)),
                       columns=['species', 'frequency'])
    f_n_out1 = os.path.join(
        file_dir, "output", "species_cycle_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])

    # load spe and reaction info
    spe_ind_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

    # convert species reaction index to real species and reactions
    d_f['species'] = d_f['species'].apply(
        lambda x: psri.pathname_to_real_spe_reaction(
            spe_ind_name_dict, new_ind_reaction_dict, x)
        .strip())
    f_n_out2 = os.path.join(file_dir, "output", "species_cycle_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])


if __name__ == "__main__":
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # species_count(FILE_DIR)
    reaction_count(FILE_DIR)
    # species_cycle(FILE_DIR)
