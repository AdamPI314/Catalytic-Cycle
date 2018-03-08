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
from math import log10
from decimal import Decimal
import networkx as nx
import parse_spe_reaction_info as psri
import interpolation
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
        if spe in dict_s:
            return dict_s[str(spe)]
        else:
            return spe
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


def rescale_array_v2(arr, min_t1=0.0, max_t1=1.0, min_t2=1.0, max_t2=10.0, threshold=None):
    """
    rescale array, to automatically select the chattering groups
    """
    arr2 = copy.deepcopy(arr)
    for idx, val in enumerate(arr2):
        if val == 0.0:
            arr2[idx] = -1000.0
        else:
            d = Decimal(val)
            # arr2[idx] = log10(val)
            arr2[idx] = float(d.log10())

    min_val = min(arr2)
    max_val = max(arr2)
    if max_val == min_val:
        for idx, val in enumerate(arr2):
            arr2[idx] = (max_t1 + min_t1) / 2.0
        return arr2
    else:
        ratio1 = (max_t1 - min_t1) / (max_val - min_val)
        ratio2 = (max_t2 - min_t2) / (max_val - min_val)
        for idx, val in enumerate(arr2):
            if val <= threshold:
                arr2[idx] = (val - min_val) * ratio1 + min_t1
            else:
                arr2[idx] = (val - min_val) * ratio2 + min_t2

        return arr2


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
        data_dir, "output", prefix + "pathway_name_candidate" + suffix + ".csv")
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
        data_dir, "output", prefix + "pathway_name_candidate" + suffix + ".csv")
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


def init_directed_network_from_concentrtion_and_reaction_rate_at_a_time(data_dir, tag="M",
                                                                        tau=10.0, end_t=0.5, end_t2=None):
    """
    init directed network
    without parallel edges
    return networkx.DiGraph
    at least time-snapshot network at a time, end_t,
    a second time can be added, but don't change network structure, 
    instead add a second attribute "weight2" to edges
    """
    s_idx_2_name, _ = psri.parse_spe_info(data_dir)
    spe_alias = read_spe_alias(os.path.join(
        data_dir, "input", "spe_alias.json"))

    time_v = np.loadtxt(os.path.join(data_dir, "output",
                                     "time_dlsode_M.csv"), dtype=float, delimiter=',')
    conc_mat = np.loadtxt(os.path.join(data_dir, "output",
                                       "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")
    rxn_rates_mat = np.loadtxt(os.path.join(data_dir, "output",
                                            "reaction_rate_dlsode_" + str(tag) + ".csv"), delimiter=",")
    # the time point where reference time tau is
    # use interpolation here
    idx_array = [i for i in range(len(time_v))]
    time_axis = int(
        round(interpolation.interp1d(time_v, idx_array, tau * end_t)))
    if time_axis >= len(time_v):
        time_axis = len(time_v) - 1

    conc_v = conc_mat[time_axis, :]
    rxn_rates_v = rxn_rates_mat[time_axis, :]

    if end_t2 is not None:
        time_axis2 = int(
            round(interpolation.interp1d(time_v, idx_array, tau * end_t2)))
        if time_axis2 >= len(time_v):
            time_axis2 = len(time_v) - 1
        rxn_rates_v2 = rxn_rates_mat[time_axis2, :]

    species_set = set()
    species_pair_weight = {}
    if end_t2 is not None:
        species_pair_weight2 = {}

    # species pairs-reactions-coefficient
    s_p_r_c = psri.parse_species_pair_reaction(data_dir)
    # print(s_p_r_c)
    for s1, s2 in s_p_r_c:
        species_set.add(int(s1))
        # print(s1, s2)
        species_set.add(int(s2))
        if (int(s1), int(s2)) not in species_pair_weight:
            species_pair_weight.update({(int(s1), int(s2)): 0.0})
        if end_t2 is not None:
            if (int(s1), int(s2)) not in species_pair_weight2:
                species_pair_weight2.update({(int(s1), int(s2)): 0.0})

        for idx in s_p_r_c[(s1, s2)]:
            r_idx = int(s_p_r_c[(s1, s2)][idx]['r_idx'])
            c1 = float(s_p_r_c[(s1, s2)][idx]['c1'])
            c2 = float(s_p_r_c[(s1, s2)][idx]['c2'])
            flux = rxn_rates_v[r_idx] * c2 / c1
            species_pair_weight[(int(s1), int(s2))] += flux

            if end_t2 is not None:
                flux2 = rxn_rates_v2[r_idx] * c2 / c1
                species_pair_weight2[(int(s1), int(s2))] += flux2

    # print(species_set)
    # print(species_pair_weight)

    edge_weight_v = []
    for idx, key in enumerate(species_pair_weight):
        edge_weight_v.append(float(species_pair_weight[key]))
    if end_t2 is not None:
        edge_weight_v2 = []
        for idx, key in enumerate(species_pair_weight2):
            edge_weight_v2.append(float(species_pair_weight2[key]))

    # rescase concentrations
    conc_v = rescale_array(conc_v, 10.0, 25.0)
    # edge_weight_v = rescale_array(edge_weight_v, 2.0, 25.0)
    edge_weight_v = rescale_array_v2(edge_weight_v, 1.5, 2.5, 15.0, 25.0, -9)
    if end_t2 is not None:
        # edge_weight_v2 = rescale_array(edge_weight_v2, 2.0, 25.0)
        edge_weight_v2 = rescale_array_v2(
            edge_weight_v2, 1.5, 2.5, 15.0, 25.0, -9)

    # final directed graph
    di_graph = nx.DiGraph()
    # add nodes first
    for _, val in enumerate(species_set):
        weight = float(conc_v[int(val)])
        node_name = change_spe_name(s_idx_2_name[str(val)], spe_alias, None)
        di_graph.add_node(node_name,
                          label=node_name, weight=weight)
    # add edges
    for idx, key in enumerate(species_pair_weight):
        src = key[0]
        dst = key[1]
        src_name = change_spe_name(s_idx_2_name[str(src)], spe_alias, None)
        dst_name = change_spe_name(s_idx_2_name[str(dst)], spe_alias, None)
        name = src_name + "," + dst_name
        if end_t2 is None:
            weight = float(edge_weight_v[idx])
            di_graph.add_edge(
                src_name, dst_name, name=name, weight=weight, weight2=weight)
        else:
            weight = float(edge_weight_v[idx])
            weight2 = float(edge_weight_v2[idx])
            di_graph.add_edge(
                src_name, dst_name, name=name, weight=weight, weight2=weight2)

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
    spe_alias_latex = read_spe_alias(os.path.join(
        data_dir, "input", "spe_alias_latex.json"))

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
    for idx, val in enumerate(labels):
        labels[idx] = change_spe_name(val, spe_alias_latex, None)
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
    a_x.set_frame_on(False)
    a_x.set_xticks([])  # this is needed for bbox_inches
    a_x.set_yticks([])

    if (path_idx == 1):
        a_x.set_title("P$_{" + str(path_idx) + "}$" + " = " +
                      "{:.6e}".format(float(pathprob)))
    else:
        a_x.set_title("P$_{" + str(path_idx) + "}$" + " = " +
                      "{:.2e}".format(float(pathprob)))

    # fig.tight_layout()
    # plt.subplots_adjust(left=0.01, right=0.9, top=0.9, bottom=0.01)
    fig.savefig(os.path.join(data_dir, "output", fig_name),
                bbox_inches='tight', dpi=500)
    # bbox_inches='tight', pad_inches=0, dpi=500)
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

    # TIME_AXIS, _ = tools.pathway_time_2_array_index(
    #     DATA_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #     species_path=G_S['species_path'], time=G_S['mc_t'])

    # TOP_N = 10
    # # PATH_NAME_SELECTED, PATH_PROB_SELECTED = get_top_n_pathway(DATA_DIR, top_n=TOP_N,
    # #                                                            suffix=SUFFIX, norm=True,
    # #                                                            species_path=G_S['species_path'],
    # #                                                            time_axis=TIME_AXIS,
    # #                                                            sort_by_p=True)
    # # PATH_IDX = [i for i in range(TOP_N)]

    # # PATH_IDX = [0, 1, 2, 3, 4, 6, 11, 44, 59, 66,
    # #             68, 93, 115, 138, 153, 165, 166, 245, 477]
    # PATH_NAME_ALL, PATH_PROB_ALL = get_top_n_pathway(DATA_DIR, top_n=G_S['top_n_p'],
    #                                                  suffix=SUFFIX, norm=True,
    #                                                  species_path=G_S['species_path'],
    #                                                  time_axis=TIME_AXIS,
    #                                                  sort_by_p=False)
    # PATH_IDX = [i for i in range(TOP_N)]
    # PATH_NAME_SELECTED = [PATH_NAME_ALL[I] for I in PATH_IDX]
    # PATH_PROB_SELECTED = [PATH_PROB_ALL[I] for I in PATH_IDX]

    # # filename without type appendix
    # NETWORK_FILENAME = PREFIX + "network" + SUFFIX

    # if not os.path.isfile(os.path.join(DATA_DIR, "output", NETWORK_FILENAME + ".gexf")):
    #     RN_OBJ = init_directed_network(
    #         DATA_DIR, path_idx=PATH_IDX, init_spe=G_S['init_s'],
    #         atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #         species_path=G_S['species_path'], time_axis=TIME_AXIS)
    #     network_to_gephi_input_file(
    #         RN_OBJ, DATA_DIR, NETWORK_FILENAME + ".gexf")

    # if os.path.isfile(os.path.join(DATA_DIR, "output", NETWORK_FILENAME + ".json")):
    #     for IDX_TMP, P_IDX in enumerate(PATH_IDX):
    #         PATHNAME = PATH_NAME_SELECTED[IDX_TMP]
    #         PATHPROB = PATH_PROB_SELECTED[IDX_TMP]

    #         plot_network(data_dir=DATA_DIR, fname=NETWORK_FILENAME + ".json",
    #                      pathname=PATHNAME, pathprob=PATHPROB,
    #                      path_idx=P_IDX + 1, end_t=G_S['end_t'], suffix=SUFFIX,
    #                      atom_followed=G_S["atom_f"], species_path=G_S['species_path'])

    END_T = 0.2
    END_T2 = 0.9
    RN_OBJ2 = init_directed_network_from_concentrtion_and_reaction_rate_at_a_time(DATA_DIR, tag="M",
                                                                                  tau=G_S['tau'],
                                                                                  end_t=END_T, end_t2=END_T2)
    if END_T2 is None:
        NETWORK_FILENAME = PREFIX + "network_all_species_" + \
            str(END_T) + SUFFIX
    else:
        NETWORK_FILENAME = PREFIX + "network_all_species_" + \
            str(END_T) + "_" + str(END_T2) + SUFFIX

    network_to_gephi_input_file(
        RN_OBJ2, DATA_DIR, NETWORK_FILENAME + ".gexf")
    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
