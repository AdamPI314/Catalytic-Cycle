import pandas as pd
import numpy as np
import re
import os
import sys


def parse_spe_info(filename):
    """
    parse species info from file= "os.path.join(file_dir, "output", "species_labelling.csv")"
    """
    line_content = np.genfromtxt(filename, dtype=str, delimiter='\n')

    matched_str = [re.findall("(\d+)\t-->\t([\w|\-|(|)]+)", line)[0]
                   for line in line_content]
    matched_str_reverse = [(x[1], x[0]) for x in matched_str]

    spe_ind_name_dict = dict(matched_str)
    spe_name_ind_dict = dict(matched_str_reverse)

    return spe_ind_name_dict, spe_name_ind_dict


# parse reactions
def parse_reaction_and_its_index(filename):
    """
    parse reaction info from file= "os.path.join(file_dir, "output", "reaction_labelling.csv")"
    """
    # load data
    line_content = np.genfromtxt(filename, dtype=str, delimiter='\n')
    matched_ind1_ind2_str = [re.findall("(\d+)\s+([-]?\d+)\s+([\w|+|=|>|<|(|)|\-|,]+)", line)[0] for line in line_content
                             if len(re.findall("(\d+)\s+([-]?\d+)\s+([\w|+|=|>|<|(|)|\-|,]+)", line)) != 0]
    # map the new old reaction index
    new_old_index_dict = dict()
    for i in range(len(matched_ind1_ind2_str)):
        new_old_index_dict.update(
            {matched_ind1_ind2_str[i][0]: str(matched_ind1_ind2_str[i][1])})

    # reactant arrow product
    reactant_arrow_product = [re.findall("([\w|+|(|)|\-|,]+)([=|>|<]+)([\w|+|(|)|\-|,]+)", ind1_ind2_reaction[2])[0]
                              for ind1_ind2_reaction in matched_ind1_ind2_str]
    reactant = [x[0] for x in reactant_arrow_product]
    arrow = [x[1] for x in reactant_arrow_product]
    product = [x[2] for x in reactant_arrow_product]
    # map reaction new reaction label and the exact reaction
    new_ind_reaction_dict = dict()
    for i in range(len(reactant_arrow_product)):
        if int(matched_ind1_ind2_str[i][1]) > 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: reactant[i] + '=>' + product[i]})
        elif int(matched_ind1_ind2_str[i][1]) < 0:
            new_ind_reaction_dict.update(
                {matched_ind1_ind2_str[i][0]: product[i] + '=>' + reactant[i]})

    return new_old_index_dict, new_ind_reaction_dict


def pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, pathway_name):
    """
    converted path to their real species name and reaction format instead of index
    """
    matched_S = re.findall("S(\d+)", pathway_name)
    matched_R = re.findall("R(\d+)", pathway_name)
    # print matched_S
    # print matched_R
    # always starts from species
    str_t = '[' + spe_ind_name_dict[matched_S[0]] + '] '
    for i in range(len(matched_R)):
        str_t += new_ind_reaction_dict[matched_R[i]]
        str_t += "-->"
        str_t += '[' + spe_ind_name_dict[matched_S[i + 1]] + '] '

    return str_t


def read_pathname_convert_2_real_spe_reaction(filename_spe, filename_reaction, filename_p, filename_p_out, topN=50):
    """
    read species and reaction info, convert path info into real species and reaction instead of index and write to file
    """
    # load path data
    path_data = pd.read_csv(filename_p, names=['path', 'prob'])
    total_prob = sum(path_data['prob'])
    # map will return a new array, will not change the value of pandas frame in situ
    # map(lambda x:x/total_prob, path_data['prob'])
    # renormalize
    path_data['prob'] /= total_prob

    # load spe and reaction info
    spe_ind_name_dict, spe_name_ind_dict = parse_spe_info(filename_spe)
    new_old_index_dict, new_ind_reaction_dict = parse_reaction_and_its_index(
        filename_reaction)

    # convert species reaction index to real species and reactions
    path_data['path'] = path_data['path'].apply(
        lambda x: pathname_to_real_spe_reaction(spe_ind_name_dict, new_ind_reaction_dict, x))

    # write to file
    path_data[0:topN].to_csv(filename_p_out, header=False,
                             index=False, sep=';', columns=['path', 'prob'])


if __name__ == '__main__':
    file_dir = os.path.abspath(os.path.realpath(
        os.path.join(sys.argv[0], os.pardir, os.pardir, os.pardir)))
    print(file_dir)
    
    # convert symbolic pathway to real pathway with real species name and real reactin name
    read_pathname_convert_2_real_spe_reaction(os.path.join(file_dir, "output", "species_labelling.csv"), os.path.join(file_dir, "output", "reaction_labelling.csv"), os.path.join(file_dir, "output", "pathway_stat.csv"), os.path.join(file_dir, "output", "pathname_prob.csv"), topN=50)
