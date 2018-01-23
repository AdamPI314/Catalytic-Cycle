"""
pattern analysis, take input from "parse pattern", do some descriptive statistics
"""

import os
import sys
from collections import OrderedDict
import numpy as np
import pandas as pd
import parse_spe_reaction_info as psri
import atom_scheme as asch
import naming
import parse_pattern
import global_settings


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
    spe_ind_name_dict, _ = psri.parse_spe_info(file_dir)
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)

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
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)
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
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)
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
    spe_ind_name_dict, _ = psri.parse_spe_info(file_dir)
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)

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
    spe_ind_name_dict, spe_name_ind_dict = psri.parse_spe_info(file_dir)
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)

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

    _, spe_name_ind_dict = psri.parse_spe_info(file_dir)
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
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(file_dir)
    # convert species reaction index to real species and reactions
    d_f['reaction'] = d_f['reaction'].apply(
        lambda x: psri.reaction_name_to_real_reaction(new_ind_reaction_dict, x)
        .strip())
    # print(d_f['reaction'])
    f_n_out2 = os.path.join(file_dir, "output", spe +
                            "_production_reaction_name.csv")
    d_f[0:top_n].to_csv(f_n_out2, header=False,
                        index=False, sep=',', columns=['reaction', 'frequency'])


def parse_spe_production_along_path(file_dir, top_n=10, spe_idx=10, init_spe=62,
                                    atom_followed="C", end_t=1.0, species_path=False,
                                    axis=0, path_branching_factor=False,
                                    s_consumption=False, s_production=True):
    """
    parse species peoduction along path, note species might not explictly shown on path
    but are side products of reaction on pathway
    if path_idx is None, use top_n path
    if path_idx is not None, instead it is a list, use only selected path, the output file name
    thereafter ends with "selected_path"
    """
    id_tmp = ""
    if spe_idx is None or spe_idx is []:
        return
    elif isinstance(spe_idx, int):
        id_tmp = str(spe_idx)
        spe_idx = [spe_idx]
    else:
        for x_t in spe_idx:
            if id_tmp == "":
                id_tmp = str(x_t)
            else:
                id_tmp += "_" + str(x_t)

    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    f_n_path_name = os.path.join(
        file_dir, "output", prefix + "pathway_name_selected" + suffix + ".csv")
    pathname_data = np.genfromtxt(
        f_n_path_name, dtype=str, max_rows=top_n + 1)

    # in case of two dimensional pathway name
    if len(np.shape(pathname_data)) == 2:
        pathname_data = pathname_data[:, axis]

    net_reactant = psri.parse_reaction_net_reactant(file_dir)
    net_product = psri.parse_reaction_net_product(file_dir)
    s_p_r_c = psri.parse_species_pair_reaction(file_dir)

    if path_branching_factor is True:
        atom_scheme = asch.get_atom_scheme(file_dir)
        s_idx_name, _ = psri.parse_spe_info(file_dir)

    s_p_c = []
    for _, p_n in enumerate(pathname_data):
        spe_consumption_count = 0
        spe_production_count = 0
        for s_i in spe_idx:
            if s_consumption is True:
                spe_consumption_count += parse_pattern.parse_species_along_path_using_reaction(
                    p_n, net_reactant, s_i, s_p_r_c)
            if s_production is True:
                spe_production_count += parse_pattern.parse_species_along_path_using_reaction(
                    p_n, net_product, s_i, s_p_r_c)

        path_branching_number = 1
        if path_branching_factor is True:
            path_branching_number = parse_pattern.calculate_path_branching_number(
                pathname=p_n, net_reactant=net_reactant, net_product=net_product,
                s_idx_name=s_idx_name, atom_scheme=atom_scheme, atom_followed=atom_followed)

        s_p_c.append((spe_production_count - spe_consumption_count)
                     * path_branching_number)

    if id_tmp != "":
        suffix += "_" + id_tmp
    f_n_spe_production_count = os.path.join(
        file_dir, "output", prefix + "pathway_species_production_count" + suffix + ".csv")

    np.savetxt(f_n_spe_production_count, s_p_c, fmt='%d')


def calculate_Merchant_alpha_value(file_dir, init_spe=10, atom_followed="C", end_t=1.0, species_path=False,
                                   s_idx=10, r_idx=736):
    """
    calculate Merchat alpha value at time point as in time.csv, not at time zero
    """
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    spe_conc_mat = np.loadtxt(os.path.join(file_dir, "output",
                                           "concentration_dlsode_M.csv"), dtype=float, delimiter=',')
    spe_k_mat = np.loadtxt(os.path.join(file_dir, "output",
                                        "drc_dlsode_M.csv"), dtype=float, delimiter=',')

    reaction_rate_mat = np.loadtxt(os.path.join(file_dir, "output",
                                                "reaction_rate_dlsode_M.csv"), dtype=float, delimiter=',')

    n_points = np.shape(spe_conc_mat)[0]
    spe_total_sink_rate_vec = np.zeros(n_points)
    merchant_alpha_v = np.zeros(n_points)

    for i in range(1, n_points):
        spe_total_sink_rate_vec[i] = spe_conc_mat[i,
                                                  s_idx] * spe_k_mat[i, s_idx]
        if spe_total_sink_rate_vec[i] > 0:
            merchant_alpha_v[i] = reaction_rate_mat[i,
                                                    r_idx] / spe_total_sink_rate_vec[i]
    # set the first value at time zero to be the same as value at time one
    merchant_alpha_v[0] = merchant_alpha_v[1]

    merchant_alpha_fn = os.path.join(
        file_dir, "output", prefix + "Merchant_alpha_" + "S" + str(s_idx) + "_R" + str(r_idx) + suffix + ".csv")
    np.savetxt(merchant_alpha_fn, merchant_alpha_v, fmt='%.15e')

    return


def parse_pathway_contains_species(file_dir, s_idx_l=None, init_spe=60,
                                   atom_followed="C", end_t=1.0, species_path=False):
    """
    parse pathway contains only species from list
    """
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    f_n_path_name = os.path.join(
        file_dir, "output", prefix + "pathway_name_selected" + suffix + ".csv")
    p_name = np.genfromtxt(f_n_path_name, dtype=str, delimiter=',')
    print(p_name)


if __name__ == "__main__":
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    G_S = global_settings.get_setting(FILE_DIR)

    # species_count(FILE_DIR)
    # reaction_count(FILE_DIR)
    # initiation_reaction_count(FILE_DIR)
    # species_cycle(FILE_DIR)
    # print(parse_species_production_path("S114R15S9R15S9", 'S9'))
    # species_production_path(FILE_DIR, spe='OH', top_n=50)
    # print(parse_species_production_reaction("S114R15S9R47S9", 'S9'))
    # species_production_reaction(FILE_DIR, spe='OH', top_n=50)
    # SPE_LIST = [14, 59, 17, 44, 38, 86,  69, 15, 82]
    # for es in SPE_LIST:
    #     path_length_statistics(
    #         FILE_DIR, init_spe=62, atom_followed="C", end_t=0.9, end_spe=es)

    # parse_spe_production_along_path(FILE_DIR, top_n=G_S['top_n_p'], spe_idx=[10],
    #                                 init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                                 end_t=G_S['end_t'], species_path=G_S['species_path'],
    #                                 axis=0, path_branching_factor=False,
    #                                 s_consumption=False, s_production=True)

    # calculate_Merchant_alpha_value(FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                                end_t=G_S['end_t'], species_path=G_S['species_path'],
    #                                s_idx=10, r_idx=736)

    parse_pathway_contains_species(FILE_DIR, s_idx_l=[60],
                                   init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
                                   end_t=G_S['end_t'], species_path=G_S['species_path'])
