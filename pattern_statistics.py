"""
pattern analysis, take input from "parse pattern", do some descriptive statistics
"""

import os
import sys
from collections import OrderedDict
import numpy as np
import pandas as pd
import parse_spe_reaction_info as psri
import naming
import parse_pattern


def path_prob_terminating_with_spe(file_dir, init_spe=62, atom_followed="C", end_t=1.0, end_spe=None, species_path=False):
    """
    get pathway and their pathway probability, path ending with spe
    """
    suffix = naming.get_suffix(file_dir=file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    if species_path is False:
        prefix = ""
    else:
        prefix = "species_"

    f_n_n = os.path.join(file_dir, "output",
                         prefix + "pathway_name_selected" + suffix + ".csv")
    f_n_p = os.path.join(file_dir, "output", prefix +
                         "pathway_prob" + suffix + ".csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    d_f = pd.DataFrame(np.transpose([pathway_name, pathway_prob]), columns=[
        'pathway', 'frequency'])

    if end_spe is None or end_spe is "ALL":
        return d_f
    else:
        d_f = d_f.loc[lambda x: x['pathway'].str.endswith("S" + str(end_spe))]
        d_f.reset_index(drop=True, inplace=True)
        # print(d_f.head())

        return d_f


def path_length_statistics(file_dir, init_spe=62, atom_followed="C", end_t=1.0, end_spe=None):
    """
    path length statistics
    """
    d_f = path_prob_terminating_with_spe(
        file_dir, init_spe, atom_followed, end_t, end_spe)

    count_map = OrderedDict()
    for _, val in enumerate(d_f['pathway'][0:20]):
        count = int(parse_pattern.parse_path_length(val))
        if count in count_map:
            count_map[count] += 1
        else:
            count_map[count] = 1
    mat = []
    count_map = OrderedDict(sorted(count_map.items()))
    for key, value in count_map.items():
        mat.append([int(key), int(value)])
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)

    if end_spe is not None:
        suffix += "_S" + str(end_spe)
    out_f_n = os.path.join(file_dir, "output", "path_length" + suffix + ".csv")
    np.savetxt(out_f_n, mat, fmt='%d', delimiter=',',
               newline='\n', header='', footer='', comments='# ')


def species_count(file_dir, top_n=50, norm=False):
    """
    species occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_species(p_n)
        for key, value in map_tmp.items():
            if key not in spe_map:
                spe_map[key] = value * p_p
            else:
                spe_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(spe_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'species', 'frequency'])

    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
    f_n_out1 = os.path.join(file_dir, "output", "species_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])

    # load spe and reaction info
    spe_ind_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "input", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))

    # convert species reaction index to real species and reactions
    d_f['species'] = d_f['species'].apply(
        lambda x: psri.pathname_to_real_spe_reaction(
            spe_ind_name_dict, new_ind_reaction_dict, x)
        .strip())
    f_n_out2 = os.path.join(file_dir, "output", "species_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])


def reaction_count(file_dir, top_n=50, norm=False):
    """
    reaction occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    reaction_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_reaction(p_n)
        for key, value in map_tmp.items():
            if key not in reaction_map:
                reaction_map[key] = value * p_p
            else:
                reaction_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(reaction_map.items(), key=lambda x: x[1], reverse=True)),
                       columns=['reaction', 'frequency'])
    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
    f_n_out1 = os.path.join(file_dir, "output", "reaction_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))
    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(file_dir, "output", "reaction_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


def initiation_reaction_count(file_dir, top_n=50, norm=False):
    """
    initiation reaction occurence in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_initiation_reaction(p_n)
        for key, value in map_tmp.items():
            if key not in spe_map:
                spe_map[key] = value * p_p
            else:
                spe_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(spe_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'reaction', 'frequency'])
    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
    f_n_out1 = os.path.join(
        file_dir, "output", "initiation_reaction_count_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))
    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(
        file_dir, "output", "initiation_reaction_count_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


def species_cycle(file_dir, top_n=50, norm=False):
    """
    species cycle in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    spe_cycle_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_species_cycle(p_n)
        for key, value in map_tmp.items():
            if key not in spe_cycle_map:
                spe_cycle_map[key] = value * p_p
            else:
                spe_cycle_map[key] += value * p_p

    d_f = pd.DataFrame(list(sorted(spe_cycle_map.items(), key=lambda x: x[1], reverse=True)),
                       columns=['species', 'frequency'])
    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
    f_n_out1 = os.path.join(
        file_dir, "output", "species_cycle_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])

    # load spe and reaction info
    spe_ind_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "input", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))

    # convert species reaction index to real species and reactions
    d_f['species'] = d_f['species'].apply(
        lambda x: psri.pathname_to_real_spe_reaction(
            spe_ind_name_dict, new_ind_reaction_dict, x)
        .strip())
    f_n_out2 = os.path.join(file_dir, "output", "species_cycle_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['species', 'frequency'])


def species_production_path(file_dir, spe='OH', top_n=50, norm=False):
    """
    species production in a path multiply by pathway probability, count pathway
    or sub-pathway ends with a species
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    # load spe and reaction info
    spe_ind_name_dict, spe_name_ind_dict = psri.parse_spe_info(os.path.join(
        file_dir, "input", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))

    species_production_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_species_production_path(
            p_n, 'S' + spe_name_ind_dict[spe])
        for key, value in map_tmp.items():
            if key not in species_production_map:
                species_production_map[key] = value * p_p
            else:
                species_production_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(species_production_map.items(), key=lambda x: x[1], reverse=True)),
                       columns=['species', 'frequency'])
    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
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


def species_production_reaction(file_dir, spe='OH', top_n=50, norm=False):
    """
    species production reaction in a path multiply by pathway probability
    """
    print(file_dir)
    f_n_n = os.path.join(file_dir, "output", "pathway_name_selected.csv")
    f_n_p = os.path.join(file_dir, "output", "pathway_prob.csv")

    pathway_name = np.genfromtxt(f_n_n, dtype=str, delimiter='\n')
    pathway_prob = np.genfromtxt(f_n_p, dtype=float, delimiter='\n')

    _, spe_name_ind_dict = psri.parse_spe_info(os.path.join(
        file_dir, "input", "species_labelling.csv"))
    reaction_map = dict()
    for _, (p_n, p_p) in enumerate(zip(pathway_name, pathway_prob)):
        map_tmp = parse_pattern.parse_species_production_reaction(
            p_n, 'S' + spe_name_ind_dict[spe])
        for key, value in map_tmp.items():
            if key not in reaction_map:
                reaction_map[key] = value * p_p
            else:
                reaction_map[key] += value * p_p
    d_f = pd.DataFrame(list(sorted(reaction_map.items(), key=lambda x: x[1], reverse=True)), columns=[
        'reaction', 'frequency'])
    if norm is True:
        total = sum(d_f['frequency'])
        d_f['frequency'] /= total
    f_n_out1 = os.path.join(file_dir, "output", spe +
                            "_production_reaction_index.csv")
    d_f[0:top_n].to_csv(f_n_out1, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])

    # load reaction info
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "input", "reaction_labelling.csv"))
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
    # species_production_path(FILE_DIR, spe='OH', top_n=50)
    # print(parse_species_production_reaction("S114R15S9R47S9", 'S9'))
    # species_production_reaction(FILE_DIR, spe='OH', top_n=50)
    SPE_LIST = [14, 59, 17, 44, 38, 86,  69, 15, 82]
    for es in SPE_LIST:
        path_length_statistics(
            FILE_DIR, init_spe=62, atom_followed="C", end_t=0.9, end_spe=es)
