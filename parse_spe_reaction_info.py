"""
parse species and reaction information
"""

import re
# import os
# import sys
import json
import pandas as pd
import numpy as np


def parse_spe_info(f_n):
    """
    parse species info from file= "os.path.join(file_dir, "output", "species_labelling.csv")"
    """
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')

    matched_str = [re.findall(r"(\d+)\t-->\t([\w|\-|(|)]+)", line)[0]
                   for line in line_content]
    matched_str_reverse = [(x[1], x[0]) for x in matched_str]

    spe_ind_name_dict = dict(matched_str)
    spe_name_ind_dict = dict(matched_str_reverse)

    return spe_ind_name_dict, spe_name_ind_dict


def read_spe_composition(f_n):
    """
    read species composition
    """
    with open(f_n, 'r') as f_h:
        data = json.load(f_h)
    return data


def parse_reaction_and_its_index(f_n):
    """
    parse reaction info from file= "os.path.join(file_dir, "output", "reaction_labelling.csv")"
    """
    # load data
    line_content = np.genfromtxt(f_n, dtype=str, delimiter='\n')
    # matched_tmp = [re.findall(r"(\d+)\s+([-]?\d+)\s+([\w|+|=|>|<|(|)|\-|\_|,]+)", line)
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
    matched_reaction = re.findall(r"R(\d+)", pathway_name)
    # only reactions
    str_t = '['
    for _, val in enumerate(matched_reaction):
        str_t += new_ind_reaction_dict[val]
    str_t += ']'
    return str_t


def pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, pathway_name):
    """
    converted path to their real species name and reaction format instead of index
    """
    matched_spe = re.findall(r"S(\d+)", pathway_name)
    matched_reaction = re.findall(r"R(\d+)", pathway_name)
    # always starts from species
    str_t = '[' + spe_ind_name_dict[matched_spe[0]] + '] '
    for idx, val in enumerate(matched_reaction):
        str_t += new_ind_reaction_dict[val]
        str_t += "-->"
        str_t += '[' + spe_ind_name_dict[matched_spe[idx + 1]] + '] '
    return str_t


def symbolic_path_2_real_path(f_n_spe, f_n_reaction, f_n_p, f_n_p_out, top_n=50, end_spe=""):
    """
    read species and reaction info,
    convert path info into real species and reaction instead of index and write to file
    """
    # load path data
    path_data = pd.read_csv(f_n_p, names=['path', 'prob'])
    total_prob = sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    # renormalize
    path_data['prob'] /= total_prob

    # filter
    if end_spe is not "" and end_spe is not "ALL":
        path_data = path_data[path_data['path'].str.endswith(end_spe)]

    # load spe and reaction info
    spe_ind_name_dict, _ = parse_spe_info(f_n_spe)
    _, new_ind_reaction_dict = parse_reaction_and_its_index(f_n_reaction)

    # convert species reaction index to real species and reactions
    path_data['path'] = path_data['path'].apply(
        lambda x: pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, x))

    # write to file
    path_data[0:top_n].to_csv(f_n_p_out, header=False,
                              index=False, sep=',', columns=['path', 'prob'])
