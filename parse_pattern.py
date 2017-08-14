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


def parse_species_production_reaction(path, spe):
    """
    parse sepcies name production, return a dictionary of species and count
    reactions ends with spe
    """
    matched_tmp = re.findall("(R\d+)" + spe, path)
    d_map = dict()
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_species_production_path(path, spe):
    """
    parse sepcies name production, return a dictionary of species and count
    pathway or sub-pathway ends with spe
    """
    idx_tmp = [(m.start(0), m.end(0)) for m in re.finditer(spe, path)]
    d_map = dict()
    # from intitial species to that species
    for _, val in enumerate(idx_tmp):
        sub_path = path[0:val[1]]
        if sub_path not in d_map:
            d_map[sub_path] = 1
        else:
            d_map[sub_path] += 1
    return d_map


def parse_reaction(path):
    """
    parse reactions, return a dictionary of reactions and their count
    """
    matched_tmp = re.findall(r"(R\d+)", path)
    d_map = dict()
    # include all reactions
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_initiation_reaction(path):
    """
    parse initiation_reactions, return a dictionary of reactions and their count
    """
    matched_tmp = re.findall(r"(R\d+)", path)
    d_map = dict()
    # include all reactions
    if len(matched_tmp) >= 1:
        for _, val in enumerate(matched_tmp[0:1]):
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

    reaction_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_reaction(p_n)
        for key, value in map_tmp.items():
            if key not in reaction_map:
                reaction_map[key] = value * p_p
            else:
                reaction_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(reaction_map.items(), key=lambda x: x[1], reverse=True)), columns=[
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
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(file_dir, "output", "reaction_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


def initiation_reaction_count(file_dir, top_n=50):
    """
    initiation reaction occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_initiation_reaction(p_n)
        for key, value in map_tmp.items():
            if key not in spe_map:
                spe_map[key] = value * p_p
            else:
                spe_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(spe_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'reaction', 'frequency'])
    f_n_out1 = os.path.join(
        file_dir, "output", "initiation_reaction_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))
    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(
        file_dir, "output", "initiation_reaction_count_name.csv")
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


def species_production_path(file_dir, spe='OH', top_n=50):
    """
    species production in a path multiply by pathway probability, count pathway
    or sub-pathway ends with a species
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    # load spe and reaction info
    spe_ind_name_dict, spe_name_ind_dict = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

    species_production_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_species_production_path(
            p_n, 'S' + spe_name_ind_dict[spe])
        for key, value in map_tmp.items():
            if key not in species_production_map:
                species_production_map[key] = value * p_p
            else:
                species_production_map[key] += value * p_p

    d_f = pd.DataFrame(list(sorted(species_production_map.items(), key=lambda x: x[1], reverse=True)),
                       columns=['species', 'frequency'])
    f_n_out1 = os.path.join(
        file_dir, "output", spe + "_production_path_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])

    # convert species reaction index to real species and reactions
    d_f['species'] = d_f['species'].apply(
        lambda x: psri.pathname_to_real_spe_reaction(
            spe_ind_name_dict, new_ind_reaction_dict, x)
        .strip())
    f_n_out2 = os.path.join(file_dir, "output", spe +
                            "_production_path_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])


def species_production_reaction(file_dir, spe='OH', top_n=50):
    """
    species production reaction in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    _, spe_name_ind_dict = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    reaction_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_species_production_reaction(
            p_n, 'S' + spe_name_ind_dict[spe])
        for key, value in map_tmp.items():
            if key not in reaction_map:
                reaction_map[key] = value * p_p
            else:
                reaction_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(reaction_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'reaction', 'frequency'])
    f_n_out1 = os.path.join(file_dir, "output", spe +
                            "_production_reaction_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))
    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(file_dir, "output", spe +
                            "_production_reaction_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


if __name__ == "__main__":
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # species_count(FILE_DIR)
    # reaction_count(FILE_DIR)
    # initiation_reaction_count(FILE_DIR)
    # species_cycle(FILE_DIR)
    # print(parse_species_production_path("S114R15S9R15S9", 'S9'))
    species_production_path(FILE_DIR, spe='OH', top_n=50)
    # print(parse_species_production_reaction("S114R15S9R47S9", 'S9'))
    species_production_reaction(FILE_DIR, spe='OH', top_n=50)
