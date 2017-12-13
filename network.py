"""
network related, generate network, save to json, save to gephi compatible file
"""

import os
import time
import sys
import re
import json
import numpy as np
import pandas as pd
import networkx as nx
import parse_spe_reaction_info as psri
import global_settings
from naming import get_suffix
import tools

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt


def read_spe_alias(filename=None):
    """
    read species alias
    """
    if filename is None:
        return None
    with open(filename, 'r') as f_h:
        spe_alias = json.load(f_h)
    return spe_alias


def update_species_idx_name_dict(dict_s, spe_alias=None):
    """
    update species idx name using alias
    """
    if spe_alias is None:
        return
    for _, val in enumerate(dict_s):
        if dict_s[val] in spe_alias.keys():
            dict_s[val] = spe_alias[dict_s[val]]
    return dict_s


def back_2_old_name(filename):
    """
    new name back to old name
    """
    spe_alias = read_spe_alias(filename)
    new_2_old = dict()
    for _, val in enumerate(spe_alias):
        new_2_old[spe_alias[val]] = val

    return new_2_old


def change_spe_name(spe, dict_s=None):
    """
    change species name, for example from S9 to H2O2
    """
    if dict_s is None:
        return spe
    return dict_s[str(spe)]


def change_rxn_name(rxn, dict_r=None):
    """
    change rxn name, for example from R9 to H2O2->OH+OH
    """
    if dict_r is None:
        return rxn
    return dict_r[str(rxn)]


def rescale_array(arr, min_t=0.0, max_t=1.0):
    """
    rescale array
    """
    min_val = min(arr)
    max_val = max(arr)
    if max_val == min_val:
        for idx, val in enumerate(arr):
            arr[idx] = (max_t + min_t) / 2.0
        return arr
    else:
        ratio = (max_t - min_t) / (max_val - min_val)
        for idx, val in enumerate(arr):
            arr[idx] = (val - min_val) * ratio + min_t
        return arr


def get_top_n_pathway(file_dir, top_n=10, init_spe=None, atom_followed=None, tau=None, pathwayEndWith=None, norm=False):
    """
    get top n path
    """
    spe_idx_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    spe_alias = read_spe_alias(os.path.join(
        file_dir, "input", "spe_alias.json"))
    spe_idx_name_dict = update_species_idx_name_dict(
        spe_idx_name_dict, spe_alias=spe_alias)

    suffix = get_suffix(file_dir, init_spe=init_spe,
                        atom_followed=atom_followed, tau=tau, pathwayEndWith=pathwayEndWith)

    f_n_path_name = os.path.join(
        file_dir, "output", "pathway_name_selected" + suffix + ".csv")
    f_n_path_prob = os.path.join(
        file_dir, "output", "pathway_prob" + suffix + ".csv")

    p_n = np.genfromtxt(f_n_path_name, dtype=str, delimiter=',')
    p_p = np.genfromtxt(f_n_path_prob, dtype=float, delimiter=',')
    p_n = p_n[0::]
    p_p = p_p[0::]
    # set the data type seperately
    d_f_n = pd.DataFrame(p_n, columns=['name'], dtype=str)
    d_f_p = pd.DataFrame(p_p, columns=['prob'], dtype=float)
    d_f = pd.concat([d_f_n, d_f_p], axis=1)

    d_f.sort_values(by='prob', ascending=False,
                    inplace=True, na_position='last')
    d_f.reset_index(drop=True, inplace=True)

    data_tmp = list(d_f['prob'])[0:top_n]
    if norm is True:
        data_tmp /= np.sum(data_tmp)

    return list(d_f['name'])[0:top_n], data_tmp


def init_directed_network(file_dir, top_n=10, init_spe=None, atom_followed=None, tau=None, pathwayEndWith=None):
    """
    init directed network
    without parallel edges
    return networkx.DiGraph
    """
    spe_idx_name_dict, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    spe_alias = read_spe_alias(os.path.join(
        file_dir, "input", "spe_alias.json"))
    spe_idx_name_dict = update_species_idx_name_dict(
        spe_idx_name_dict, spe_alias=spe_alias)

    suffix = get_suffix(file_dir, init_spe=init_spe,
                        atom_followed=atom_followed, tau=tau, pathwayEndWith=pathwayEndWith)

    f_n_path_name = os.path.join(
        file_dir, "output", "pathway_name_selected" + suffix + ".csv")
    f_n_path_prob = os.path.join(
        file_dir, "output", "pathway_prob" + suffix + ".csv")
    print(f_n_path_name, f_n_path_prob)
    p_n = np.genfromtxt(f_n_path_name, dtype=str, delimiter=',')
    p_p = np.genfromtxt(f_n_path_prob, dtype=float, delimiter=',')
    p_n = p_n[1::]
    p_p = p_p[1::]
    # set the data type seperately
    d_f_n = pd.DataFrame(p_n, columns=['name'], dtype=str)
    d_f_p = pd.DataFrame(p_p, columns=['prob'], dtype=float)
    d_f = pd.concat([d_f_n, d_f_p], axis=1)

    d_f.sort_values(by='prob', ascending=False,
                    inplace=True, na_position='last')
    d_f.reset_index(drop=True, inplace=True)
    print(d_f.head())
    d_f = d_f.loc[0:(top_n + 1)]

    # temporary directed graph
    d_g_tmp = nx.DiGraph()

    # record all nodes
    nodes = set()
    for _, val in d_f.iterrows():
        matched_spe = re.findall(r"S(\d+)", val['name'])
        for _, spe in enumerate(matched_spe):
            nodes.add(change_spe_name(spe, spe_idx_name_dict))

    for _, val in enumerate(nodes):
        d_g_tmp.add_node(val, weight=0.0, label=str(val))

    for _, val in d_f.iterrows():
        prob = float(val['prob'])

        # get rid of R-1000003S90, don't need it here
        print(val['name'])
        path_name_tmp = re.sub(r"R-\d+S\d+", r'', val['name'])
        print(path_name_tmp)

        matched_spe = re.findall(r"S(\d+)", path_name_tmp)
        matched_reaction = re.findall(r"R(\d+)", path_name_tmp)
        for idx, spe in enumerate(matched_spe):
            d_g_tmp.node[change_spe_name(
                spe, spe_idx_name_dict)]['weight'] += 1.0 * prob
            if idx > 0:
                src = change_spe_name(matched_spe[idx - 1], spe_idx_name_dict)
                dest = change_spe_name(spe, spe_idx_name_dict)
                rxn = change_rxn_name(matched_reaction[idx - 1])
                if d_g_tmp.has_edge(src, dest):
                    d_g_tmp[src][dest]['weight'] += 1.0 * prob
                    d_g_tmp[src][dest]['reactions'].add(rxn)
                else:
                    d_g_tmp.add_edge(src, dest,
                                     reactions=set([rxn]), weight=1.0 * prob)

    # update directed graph, for example,
    # 1. reactions is originally a set, combine to get a string of reactions
    # 2. smooth and re-normalize node weight
    # 3. re-normalize edge weight
    node_weight = []
    for _, val in enumerate(d_g_tmp.nodes()):
        node_weight.append(d_g_tmp.node[val]['weight'])
    edge_weight = []
    for _, val in enumerate(d_g_tmp.edges()):
        edge_weight.append(d_g_tmp[val[0]][val[1]]['weight'])
    node_weight = rescale_array(node_weight, 1.0, 5.0)
    edge_weight = rescale_array(edge_weight, 3.0, 15.0)
    # final directed graph
    di_graph = nx.DiGraph()
    for idx, val in enumerate(d_g_tmp.nodes()):
        di_graph.add_node(val, weight=node_weight[idx])
    for idx, val in enumerate(d_g_tmp.edges()):
        src = val[0]
        dest = val[1]

        rxn_set = d_g_tmp[src][dest]['reactions']
        rxn_set = sorted(rxn_set, key=lambda x: int(x), reverse=False)
        name = ",".join(x for x in rxn_set)

        weight = edge_weight[idx]
        di_graph.add_edge(src, dest, name=name, weight=weight)

    return di_graph


def network_to_gephi_input_file(networkx_obj, file_dir, fname="network.gexf"):
    """
    write to gephi input file
    """
    nx.write_gexf(networkx_obj, os.path.join(file_dir, "output", fname))
    return


def get_names_coordinates(file_dir, fname=""):
    """
    read nodes names and corresponding coordinates
    """
    with open(os.path.join(file_dir, "output", fname), 'r') as f_h:
        data = json.load(f_h)

    # naming thing
    new_2_old = back_2_old_name(os.path.join(
        file_dir, "input", "spe_alias.json"))

    # we only need names and coordinates
    n_coordinate = dict()
    for _, node in enumerate(data['nodes']):
        if node['label'] in new_2_old:
            n_coordinate[new_2_old[node['label']]] = (node['x'], node['y'])
        else:
            n_coordinate[node['label']] = (node['x'], node['y'])

    # print(n_coordinate)
    return n_coordinate


def plot_network(file_dir, fname="", pathname="", pathprob=1.0, flag="", tau=1.0):
    """
    plot network manually
    """
    print(fname)
    n_coordinate = get_names_coordinates(file_dir, fname)

    # figure name
    if flag is "":
        fig_name = "network_path" + ".jpg"
    else:
        fig_name = "network_path_" + str(flag) + ".jpg"

    # specify label for lines
    labels = []
    x = []
    y = []
    name_idx_dict = dict()
    for i_tmp, val in enumerate(n_coordinate):
        labels.append(val)
        name_idx_dict[val] = i_tmp
        x.append(float(n_coordinate[val][0]))
        y.append(float(n_coordinate[val][1]))

    # read in species index name
    spe_idx_name_dict, spe_name_idx_dict = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

    # modify labels
    spe_union_find_group = global_settings.get_union_find_group(FILE_DIR)
    for idx, val in enumerate(labels):
        spe_i = spe_name_idx_dict[val]
        if spe_i in spe_union_find_group:
            labels[idx] = ",".join([str(spe_idx_name_dict[str(x)])
                                    for x in spe_union_find_group[spe_i]])
    print(labels)

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    # background
    a_x.scatter(x, y,
                color='b', marker="o", alpha=0.3)
    for i, _ in enumerate(x):
        t_h = a_x.annotate(labels[i], (x[i], y[i]))
        t_h.set_alpha(0.15)

    # get rid of R-1000003S90, don't need it here
    pathname = re.sub(r"R-\d+S\d+", r'', pathname)
    # parse pathway
    matched_spe = re.findall(r"S(\d+)", pathname)
    matched_reaction = re.findall(r"R(\d+)", pathname)
    print(matched_reaction)
    node_list = [name_idx_dict[spe_idx_name_dict[str(x)]] for x in matched_spe]
    print(node_list)
    for idx, curr_idx in enumerate(node_list):
        if idx >= 1:
            pre_idx = node_list[idx - 1]
            a_h = a_x.annotate('', xy=(x[curr_idx], y[curr_idx]), xytext=(x[pre_idx], y[pre_idx]),
                               arrowprops={'arrowstyle': '->',
                                           'lw': 4, 'color': 'red'},
                               va='center')
            a_h.set_alpha(0.9)

    # re-draw points and labels on canvas
    for _, val in enumerate(node_list):
        a_x.scatter(x[val], y[val],
                    color='r', marker="o", alpha=0.9)
        t_h = a_x.annotate(labels[val], (x[val], y[val]))
        t_h.set_alpha(0.9)

    # draw reaction along path
    for idx, curr_idx in enumerate(node_list):
        if idx >= 1:
            pre_idx = node_list[idx - 1]
            rxn_name = str(new_ind_reaction_dict[matched_reaction[idx - 1]])

            t_h = a_x.annotate(rxn_name,
                               (x[pre_idx], y[pre_idx] * 0.5 + y[curr_idx] * 0.5), color='g', size=8.0)
            t_h.set_alpha(0.5)

    a_x.set_xlim([np.min(x) - 0.01 * (np.max(x) - np.min(x)),
                  np.max(x) + 0.25 * (np.max(x) - np.min(x))])
    # a_x.grid('on')
    a_x.axis('off')
    a_x.set_title(flag + "(" + str(tau) + "$\\tau$" + ")" +
                  " = " + "{:.2e}".format(float(pathprob)))

    # fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()

    return


if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    G_S = global_settings.get_setting(FILE_DIR)

    ATOM_FOLLOWED = G_S['atom_f']
    # PREFIX = "C3H8"
    PREFIX = "S" + str(G_S['init_s'])

    # RN_OBJ = init_directed_network(
    #     FILE_DIR, top_n=G_S['top_n_p_gephi'], init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], tau=G_S['tau'], pathwayEndWith=None)
    # network_to_gephi_input_file(
    #     RN_OBJ, FILE_DIR, PREFIX + "_" + G_S['atom_f'] + "_network_" + str(G_S['top_n_p_gephi']) + "_" + str(G_S['tau']) + ".gexf")

    PATH_NAME_TOP_N, PATH_PROB_TOP_N = get_top_n_pathway(
        FILE_DIR, top_n=50, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
        tau=G_S['tau'], pathwayEndWith=None, norm=True)
    for idx, pathname in enumerate(PATH_NAME_TOP_N):
        plot_network(file_dir=FILE_DIR, fname=PREFIX + "_" +
                     G_S['atom_f'] + "_network_" +
                     str(G_S['top_n_p_gephi']) + "_" +
                     str(G_S['tau']) + ".json",
                     pathname=pathname, pathprob=PATH_PROB_TOP_N[idx], flag="P" + str(idx + 1), tau=G_S['tau'])

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
