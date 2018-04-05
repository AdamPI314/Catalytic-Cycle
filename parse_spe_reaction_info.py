"""
parse species and reaction information
"""

import re
import os
import sys
import json
import pandas as pd
import numpy as np
import read_write_configuration as rwc
# import time


def parse_spe_info(data_dir):
    """
    parse species info from file= "os.path.join(data_dir, "output", "species_labelling.csv")"
    """
    f_n = os.path.join(data_dir, "input", "species_labelling.csv")
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')

    matched_str = [re.findall(r"(\d+)\t-->\t([\w|\-|(|)]+)", line)[0]
                   for line in line_content]
    matched_str_reverse = [(x[1], x[0]) for x in matched_str]

    spe_ind_name_dict = dict(matched_str)
    spe_name_ind_dict = dict(matched_str_reverse)

    return spe_ind_name_dict, spe_name_ind_dict


def parse_species_pair_reaction(data_dir):
    """
    parse species pairs and associated reactions, coefficients
    """
    f_n = os.path.join(data_dir, "input", "species_pairs_reactions_coefs.json")

    s_p_r_c_old = rwc.read_configuration(f_n)
    s_p_r_c_pair_counter = {}
    s_p_r_c_new = {}
    for idx in s_p_r_c_old:
        src = s_p_r_c_old[idx]['from']
        dst = s_p_r_c_old[idx]['to']
        r_idx = s_p_r_c_old[idx]['r_idx']
        c1 = s_p_r_c_old[idx]['c1']
        c2 = s_p_r_c_old[idx]['c2']
        r_name = s_p_r_c_old[idx]['r_name']

        base_d = {'r_idx': r_idx, 'c1': c1, 'c2': c2, 'r_name': r_name}
        if (src, dst) not in s_p_r_c_new:
            s_p_r_c_pair_counter[(src, dst)] = 1
            s_p_r_c_new[(src, dst)] = {'0': base_d}
        else:
            counter = s_p_r_c_pair_counter[(src, dst)]
            s_p_r_c_new[(src, dst)].update({str(counter): base_d})
            s_p_r_c_pair_counter[(src, dst)] += 1

    # for key in s_p_r_c_new:
        # print(key)
        # print(s_p_r_c_new[key])
    return s_p_r_c_new


def get_reactions_from_si_2_sj(data_dir, si, sj):
    """
    return a list of reactions that transform species si to sj
    """
    s_p_r_c = parse_species_pair_reaction(data_dir)
    r_list = []
    if (str(si), str(sj)) in s_p_r_c:
        for _, val in s_p_r_c[(str(si), str(sj))].items():
            r_list.append(int(val['r_idx']))
    print(r_list)
    return r_list


def read_spe_composition(f_n):
    """
    read species composition
    """
    with open(f_n, 'r') as f_h:
        data = json.load(f_h)
    return data


def parse_reaction_and_its_index(data_dir):
    """
    parse reaction info from file= "os.path.join(data_dir, "input", "reaction_labelling.csv")"
    """
    f_n = os.path.join(data_dir, "input", "reaction_labelling.csv")
    # load data
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')
    matched_tmp = [re.findall(r"([\d]+)\s+([\-\d]+)\s+([\w\(\)\-\_,\+]+\<?={1}\>?[\w\(\)\-\_,\+]+)", line)
                   for line in line_content]
    matched_ind1_ind2_str = [x[0] for x in matched_tmp if len(x) != 0]
    # map the new old reaction index
    new_old_index_dict = dict()
    for _, val in enumerate(matched_ind1_ind2_str):
        new_old_index_dict.update(
            {val[0]: str(val[1])})

    # reactant arrow product
    reactant_product = [re.findall(r"([\w|+|(|)|\-|\_|,]+)[=|>|<]+([\w|+|(|)|\-|\_|,]+)",
                                   ind1_ind2_reaction[2])[0]
                        for ind1_ind2_reaction in matched_ind1_ind2_str]
    reactant = [x[0] for x in reactant_product]
    product = [x[1] for x in reactant_product]
    # map reaction new reaction label and the exact reaction
    new_ind_reaction_dict = dict()
    for i in range(len(reactant_product)):
        if int(matched_ind1_ind2_str[i][1]) > 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: reactant[i] + '=>' + product[i]})
        elif int(matched_ind1_ind2_str[i][1]) < 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: product[i] + '=>' + reactant[i]})
    return new_old_index_dict, new_ind_reaction_dict


def reaction_name_to_real_reaction(new_ind_reaction_dict, pathway_name):
    """
    converted reaction name to their reaction format instead of index
    """
    matched_reaction = re.findall(r"R([-]?\d+)", pathway_name)
    # only reactions
    str_t = '['
    for _, val in enumerate(matched_reaction):
        if '-' not in val:
            str_t += new_ind_reaction_dict[val]
        else:
            str_t += "<-chattering->"
    str_t += ']'
    return str_t


def pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, pathway_name):
    """
    converted path to their real species name and reaction format instead of index
    """
    # always starts from species
    str_t = ""
    matched_s_r = re.findall(r"S\d+(?:R[-]?\d+)?", pathway_name)
    for idx, val in enumerate(matched_s_r):
        m_s = re.findall(r"S(\d+)", val)
        m_r = re.findall(r"R([-]?\d+)", val)

        m_s_idx = m_s[0]
        str_t += '[' + spe_ind_name_dict[m_s_idx] + ']'
        if (len(m_r) == 0):
            if idx != len(matched_s_r) - 1:
                str_t += "-->"
        elif len(m_r) > 0:
            m_r_idx = m_r[0]
            if '-' not in m_r_idx:
                str_t += new_ind_reaction_dict[m_r_idx]
                str_t += "-->"
            else:
                str_t += "<-chattering->"

    return str_t


def symbolic_path_2_real_path(data_dir, f_n_p, f_n_p_out, top_n=50, end_s_idx=None, max_rows=5000):
    """
    read species and reaction info,
    convert path info into real species and reaction instead of index and write to file
    """

    # load path data
    if end_s_idx is None or end_s_idx is []:
        path_data = pd.read_csv(f_n_p, names=['path', 'prob'], nrows=top_n + 1)
    elif end_s_idx is not None and end_s_idx is not []:
        n_spe = len(end_s_idx)
        path_data = pd.read_csv(
            f_n_p, names=['path', 'prob'], nrows=top_n * n_spe + max_rows)

    total_prob = sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    # renormalize
    path_data['prob'] /= total_prob

    # filter
    if end_s_idx is not None and end_s_idx is not []:
        end_spe_str = ['S' + str(x) for x in end_s_idx]
        end_spe_tuple = tuple(end_spe_str)
        path_data = path_data[path_data['path'].str.endswith(end_spe_tuple)]

    # load spe and reaction info
    spe_ind_name_dict, _ = parse_spe_info(data_dir)
    _, new_ind_reaction_dict = parse_reaction_and_its_index(data_dir)

    # convert species reaction index to real species and reactions
    path_data['path'] = path_data['path'].apply(
        lambda x: pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, x))

    # write to file
    path_data[0:top_n].to_csv(f_n_p_out, header=False,
                              index=False, sep=',', columns=['path', 'prob'])


def symbolic_path_2_real_path_pff(data_dir, fn):
    """
    convert symbolic pathway to real pathway with real species name and real reaction name
    path from file
    """
    out_fn = fn[0:-4] + "_real_path" + ".csv"

    symbolic_path_2_real_path(
        data_dir,
        os.path.join(
            data_dir, "output", fn),
        os.path.join(
            data_dir, "output", out_fn),
        10000000, None)


def parse_reaction_net_reactant(data_dir):
    """
    return a dict of "species": number based on reaction reactant
    """
    f_n = os.path.join(data_dir, "input", "reaction_information.json")

    data = rwc.read_configuration(f_n)

    net_reactant = {}
    for _, r_idx in enumerate(data):
        entry = {}
        for val1 in data[r_idx]['net_reactant']:
            entry.update({data[r_idx]['net_reactant'][val1]['species_index']:
                          data[r_idx]['net_reactant'][val1]['coefficient']})
        net_reactant.update({r_idx: entry})

    return net_reactant


def parse_reaction_net_product(data_dir):
    """
    return a dict of "species": number based on reaction product
    """
    f_n = os.path.join(data_dir, "input", "reaction_information.json")

    data = rwc.read_configuration(f_n)

    net_product = {}
    for _, r_idx in enumerate(data):
        entry = {}
        for val1 in data[r_idx]['net_product']:
            entry.update({data[r_idx]['net_product'][val1]['species_index']:
                          data[r_idx]['net_product'][val1]['coefficient']})
        net_product.update({r_idx: entry})

    return net_product


def net_source_reaction_of_species(data_dir, spe_idx):
    """
    return a list of reaction, to each reaction of this list,
    they must return at least one expected species
    """
    net_product = parse_reaction_net_product(data_dir)
    reaction_list = []
    _, rxn_i_2_n = parse_reaction_and_its_index(data_dir)
    name_list = []
    for rxn in net_product:
        if str(spe_idx) in net_product[rxn]:
            reaction_list.append(int(rxn))
            name_list.append(rxn_i_2_n[str(rxn)])

    print(reaction_list, len(reaction_list))
    print(name_list, len(name_list))
    return sorted(reaction_list)


def net_sink_reaction_of_species(data_dir, spe_idx):
    """
    return a list of reaction, to each reaction of this list,
    they are net sink reaction of a species
    """
    net_reactant = parse_reaction_net_reactant(data_dir)
    reaction_list = []
    _, rxn_i_2_n = parse_reaction_and_its_index(data_dir)
    name_list = []
    for rxn in net_reactant:
        if str(spe_idx) in net_reactant[rxn]:
            reaction_list.append(int(rxn))
            name_list.append(rxn_i_2_n[str(rxn)])

    print(reaction_list, len(reaction_list))
    print(name_list, len(name_list))
    return sorted(reaction_list)


if __name__ == '__main__':
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)

    # parse_reaction_net_product(DATA_DIR)
    # net_source_reaction_of_species(DATA_DIR, 10)
    # net_sink_reaction_of_species(DATA_DIR, 62)
    # net_sink_reaction_of_species(DATA_DIR, 10)
    # net_sink_reaction_of_species(DATA_DIR, 12)
    net_sink_reaction_of_species(DATA_DIR, 17)
    # net_sink_reaction_of_species(DATA_DIR, 16)
    # get_reactions_from_si_2_sj(DATA_DIR, 60, 78)
    # get_reactions_from_si_2_sj(DATA_DIR, 78, 60)
    # get_reactions_from_si_2_sj(DATA_DIR, 78, 87)
    # get_reactions_from_si_2_sj(DATA_DIR, 87, 78)
    # get_reactions_from_si_2_sj(DATA_DIR, 87, 90)
    # get_reactions_from_si_2_sj(DATA_DIR, 90, 87)
