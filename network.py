"""
network related, generate network, save to json, save to gephi compatible file
"""

import os
import time
import sys
import re
import json
import copy
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
    if not os.path.isfile(filename):
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
    if spe_alias is None:
        return {}
    new_2_old = dict()
    for _, val in enumerate(spe_alias):
        new_2_old[spe_alias[val]] = val

    return new_2_old


def change_spe_name(spe, dict_s=None, union_find=None):
    """
    change species name, for example from S9 to H2O2
    """
    # dict_s cannot be None, otherwise return
    if dict_s is None:
        return spe
    if union_find is None or str(spe) not in union_find:
        return dict_s[str(spe)]
    # if union find is not None and str(spe) in union_find
    # regard the first element in the set as the root element
    set_tmp = copy.deepcopy(union_find[str(spe)])
    return dict_s[str(next(iter(set_tmp)))]


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


def get_top_n_pathway(data_dir, top_n=10, suffix="", norm=False, species_path=False, time_axis=0, sort_by_p=False):
    """
    get top n path
    """
    spe_idx_name_dict, _ = psri.parse_spe_info(data_dir)
    spe_alias = read_spe_alias(os.path.join(
        data_dir, "input", "spe_alias.json"))
    spe_idx_name_dict = update_species_idx_name_dict(
        spe_idx_name_dict, spe_alias=spe_alias)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    f_n_path_name = os.path.join(
        data_dir, "output", prefix + "pathway_name_selected" + suffix + ".csv")
    f_n_path_prob = os.path.join(
        data_dir, "output", prefix + "pathway_prob" + suffix + ".csv")

    p_n = np.genfromtxt(f_n_path_name, dtype=str, delimiter=',')
    p_p = np.genfromtxt(f_n_path_prob, dtype=float, delimiter=',')

    # in case of two dimensional pathway name
    if len(np.shape(p_n)) == 2:
        p_n = p_n[:, time_axis]
    if len(np.shape(p_p)) == 2:
        p_p = p_p[:, time_axis]

    # set the data type seperately
    d_f_n = pd.DataFrame(p_n, columns=['name'], dtype=str)
    d_f_p = pd.DataFrame(p_p, columns=['prob'], dtype=float)
    d_f = pd.concat([d_f_n, d_f_p], axis=1)

    if sort_by_p is True:
        d_f.sort_values(by='prob', ascending=False,
                        inplace=True, na_position='last')
    d_f.reset_index(drop=True, inplace=True)

    data_tmp = list(d_f['prob'])[0:top_n]
    if norm is True:
        data_tmp /= np.sum(data_tmp)

    return list(d_f['name'])[0:top_n], data_tmp


def init_directed_network(data_dir, path_idx=None, init_spe=None, atom_followed="C",
                          end_t=None, species_path=False, time_axis=0):
    """
    init directed network
    without parallel edges
    return networkx.DiGraph
    """
    spe_idx_name_dict, _ = psri.parse_spe_info(data_dir)

    suffix = get_suffix(data_dir, init_spe=init_spe,
                        atom_followed=atom_followed, end_t=end_t)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    f_n_path_name = os.path.join(
        data_dir, "output", prefix + "pathway_name_selected" + suffix + ".csv")
    f_n_path_prob = os.path.join(
        data_dir, "output", prefix + "pathway_prob" + suffix + ".csv")

    print(f_n_path_name, f_n_path_prob)
    p_n = np.genfromtxt(f_n_path_name, dtype=str, delimiter=',')
    p_p = np.genfromtxt(f_n_path_prob, dtype=float, delimiter=',')

    # in case of two dimensional pathway name
    if len(np.shape(p_n)) == 2:
        p_n = p_n[:, time_axis]
    if len(np.shape(p_p)) == 2:
        p_p = p_p[:, time_axis]

    # retrieve pathway name and pathway probability before sort
    p_n = [p_n[i] for i in path_idx]
    p_p = [p_p[i] for i in path_idx]

    # set the data type seperately
    d_f_n = pd.DataFrame(p_n, columns=['name'], dtype=str)
    d_f_p = pd.DataFrame(p_p, columns=['prob'], dtype=float)
    d_f = pd.concat([d_f_n, d_f_p], axis=1)

    d_f.sort_values(by='prob', ascending=False,
                    inplace=True, na_position='last')
    d_f.reset_index(drop=True, inplace=True)
    print(d_f.head())

    # temporary directed graph
    d_g_tmp = nx.DiGraph()

    # modify labels
    spe_union_find_group = global_settings.get_union_find_group(
        DATA_DIR, atom_followed)

    # record all nodes
    nodes = set()
    for _, val in d_f.iterrows():
        matched_spe = re.findall(r"S(\d+)", val['name'])
        for _, spe in enumerate(matched_spe):
            nodes.add(change_spe_name(spe, spe_idx_name_dict,
                                      union_find=spe_union_find_group))

    for _, val in enumerate(nodes):
        d_g_tmp.add_node(val, weight=0.0, label=str(val))

    for _, val in d_f.iterrows():
        prob = float(val['prob'])

        # get rid of R-1000003S90, don't need it here
        print(val['name'])
        path_name_tmp = re.sub(r"R-\d+S\d+", r'', val['name'])
        print(path_name_tmp)

        # pathway contains both reaction and species
        if species_path is False:
            matched_spe = re.findall(r"S(\d+)", path_name_tmp)
            matched_reaction = re.findall(r"R(\d+)", path_name_tmp)
            for idx, spe in enumerate(matched_spe):
                d_g_tmp.node[change_spe_name(
                    spe, spe_idx_name_dict, union_find=spe_union_find_group)]['weight'] += 1.0 * prob
                if idx > 0:
                    src = change_spe_name(
                        matched_spe[idx - 1], spe_idx_name_dict, union_find=spe_union_find_group)
                    dest = change_spe_name(
                        spe, spe_idx_name_dict, union_find=spe_union_find_group)
                    rxn = change_rxn_name(matched_reaction[idx - 1])
                    if d_g_tmp.has_edge(src, dest):
                        d_g_tmp[src][dest]['weight'] += 1.0 * prob
                        d_g_tmp[src][dest]['reactions'].add(rxn)
                    else:
                        d_g_tmp.add_edge(src, dest,
                                         reactions=set([rxn]), weight=1.0 * prob)
        else:
            matched_spe = re.findall(r"S(\d+)", path_name_tmp)
            for idx, spe in enumerate(matched_spe):
                d_g_tmp.node[change_spe_name(
                    spe, spe_idx_name_dict, union_find=spe_union_find_group)]['weight'] += 1.0 * prob
                if idx > 0:
                    src = change_spe_name(
                        matched_spe[idx - 1], spe_idx_name_dict, union_find=spe_union_find_group)
                    dest = change_spe_name(
                        spe, spe_idx_name_dict, union_find=spe_union_find_group)
                    rxn = '-1'
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


def network_to_gephi_input_file(networkx_obj, data_dir, fname="network.gexf"):
    """
    write to gephi input file
    """
    nx.write_gexf(networkx_obj, os.path.join(data_dir, "output", fname))
    return


def get_names_coordinates(data_dir, fname=""):
    """
    read nodes names and corresponding coordinates
    """
    with open(os.path.join(data_dir, "output", fname), 'r') as f_h:
        data = json.load(f_h)

    # naming thing
    new_2_old = back_2_old_name(os.path.join(
        data_dir, "input", "spe_alias.json"))

    # we only need names and coordinates
    n_coordinate = dict()
    for _, node in enumerate(data['nodes']):
        if node['label'] in new_2_old:
            n_coordinate[new_2_old[node['label']]] = (node['x'], node['y'])
        else:
            n_coordinate[node['label']] = (node['x'], node['y'])

    # print(n_coordinate)
    return n_coordinate


def plot_network(data_dir, fname="", pathname="", pathprob=1.0, path_idx=None, end_t=1.0, suffix="", atom_followed="C", species_path=False):
    """
    plot network manually
    """
    print(fname)
    n_coordinate = get_names_coordinates(data_dir, fname)

    prefix = ""
    if species_path is True:
        prefix = "species_"

    # figure name
    if suffix is "":
        fig_name = prefix + "network_path_" + str(path_idx) + ".jpg"
    else:
        fig_name = prefix + "network_path_" + \
            str(path_idx) + str(suffix) + ".jpg"

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
    spe_idx_name_dict, spe_name_idx_dict = psri.parse_spe_info(data_dir)
    _, new_ind_reaction_dict = psri.parse_reaction_and_its_index(data_dir)

    # modify labels
    spe_union_find_group = global_settings.get_union_find_group(
        DATA_DIR, atom_followed)
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
    print(matched_spe, matched_reaction)
    node_list = [name_idx_dict[change_spe_name(
        str(x), spe_idx_name_dict, spe_union_find_group)] for x in matched_spe]
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
    if species_path is False:
        # check for duplicate transition
        idx_label_dict = {}
        for idx, curr_idx in enumerate(node_list):
            if idx >= 1:
                pre_idx = node_list[idx - 1]
                rxn_idx = matched_reaction[idx - 1]
                if tuple([pre_idx, curr_idx, rxn_idx]) in idx_label_dict:
                    idx_label_dict[tuple([pre_idx, curr_idx, rxn_idx])
                                   ] += "," + str(idx)
                else:
                    idx_label_dict[tuple(
                        [pre_idx, curr_idx, rxn_idx])] = str(idx)

        for idx, curr_idx in enumerate(node_list):
            if idx >= 1:
                pre_idx = node_list[idx - 1]
                rxn_idx = matched_reaction[idx - 1]
                rxn_name = idx_label_dict[tuple([pre_idx, curr_idx, rxn_idx])] + ": " + str(
                    new_ind_reaction_dict[matched_reaction[idx - 1]])

                if x[pre_idx] <= x[curr_idx]:
                    x_tmp = x[pre_idx]
                else:
                    x_tmp = x[curr_idx]
                y_tmp = y[pre_idx] * 0.7 + y[curr_idx] * 0.3

                t_h = a_x.annotate(rxn_name,
                                   (x_tmp, y_tmp), color='g', size=8.0)
                t_h.set_alpha(0.5)
    else:
        # build idx->label
        idx_label_dict = {}
        for idx, curr_idx in enumerate(node_list):
            if idx >= 1:
                pre_idx = node_list[idx - 1]
                if tuple([pre_idx, curr_idx]) in idx_label_dict:
                    idx_label_dict[tuple([pre_idx, curr_idx])
                                   ] += "," + str(idx)
                else:
                    idx_label_dict[tuple([pre_idx, curr_idx])] = str(idx)

        for idx, curr_idx in enumerate(node_list):
            if idx >= 1:
                pre_idx = node_list[idx - 1]
                rxn_name = idx_label_dict[tuple([pre_idx, curr_idx])]

                t_h = a_x.annotate(rxn_name,
                                   (x[pre_idx] * 0.7 + x[curr_idx] * 0.3, y[pre_idx] * 0.7 + y[curr_idx] * 0.3), color='g', size=8.0)
                t_h.set_alpha(0.5)

    a_x.set_xlim([np.min(x) - 0.01 * (np.max(x) - np.min(x)),
                  np.max(x) + 0.25 * (np.max(x) - np.min(x))])
    # a_x.grid('on')
    a_x.axis('off')
    a_x.set_title("P$_{" + str(path_idx) + "}$" + " = " +
                  "{:.2e}".format(float(pathprob)))

    # fig.tight_layout()
    fig.savefig(os.path.join(data_dir, "output", fig_name), dpi=500)
    plt.close()

    return


if __name__ == '__main__':
    INIT_TIME = time.time()

    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)

    G_S = global_settings.get_setting(DATA_DIR)

    PREFIX = ""
    if G_S['species_path'] is True:
        PREFIX = "species_"
    SUFFIX = get_suffix(DATA_DIR, init_spe=G_S['init_s'],
                        atom_followed=G_S['atom_f'], end_t=G_S['end_t'])

    TIME_AXIS, _ = tools.pathway_time_2_array_index(
        DATA_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
        species_path=G_S['species_path'], time=G_S['mc_t'])

    TOP_N = 15
    # PATH_NAME_SELECTED, PATH_PROB_SELECTED = get_top_n_pathway(DATA_DIR, top_n=TOP_N,
    #                                                            suffix=SUFFIX, norm=True,
    #                                                            species_path=G_S['species_path'],
    #                                                            time_axis=TIME_AXIS,
    #                                                            sort_by_p=True)
    # PATH_IDX = [i for i in range(TOP_N)]

    # PATH_IDX = [0, 1, 2, 3, 4, 6, 11, 44, 59, 66,
    #             68, 93, 115, 138, 153, 165, 166, 245, 477]
    PATH_NAME_ALL, PATH_PROB_ALL = get_top_n_pathway(DATA_DIR, top_n=G_S['top_n_p'],
                                                     suffix=SUFFIX, norm=True,
                                                     species_path=G_S['species_path'],
                                                     time_axis=TIME_AXIS,
                                                     sort_by_p=False)
    PATH_IDX = [i+1 for i in range(TOP_N-1)]
    PATH_NAME_SELECTED = [PATH_NAME_ALL[I] for I in PATH_IDX]
    PATH_PROB_SELECTED = [PATH_PROB_ALL[I] for I in PATH_IDX]

    # filename without type appendix
    NETWORK_FILENAME = PREFIX + "network" + SUFFIX

    if not os.path.isfile(os.path.join(DATA_DIR, "output", NETWORK_FILENAME + ".gexf")):
        RN_OBJ = init_directed_network(
            DATA_DIR, path_idx=PATH_IDX, init_spe=G_S['init_s'],
            atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
            species_path=G_S['species_path'], time_axis=TIME_AXIS)
        network_to_gephi_input_file(
            RN_OBJ, DATA_DIR, NETWORK_FILENAME + ".gexf")

    if os.path.isfile(os.path.join(DATA_DIR, "output", NETWORK_FILENAME + ".json")):
        for IDX_TMP, P_IDX in enumerate(PATH_IDX):
            PATHNAME = PATH_NAME_SELECTED[IDX_TMP]
            PATHPROB = PATH_PROB_SELECTED[IDX_TMP]

            plot_network(data_dir=DATA_DIR, fname=NETWORK_FILENAME + ".json",
                         pathname=PATHNAME, pathprob=PATHPROB,
                         path_idx=P_IDX + 1, end_t=G_S['end_t'], suffix=SUFFIX,
                         atom_followed=G_S["atom_f"], species_path=G_S['species_path'])

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
