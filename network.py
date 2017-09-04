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
    ratio = (max_t - min_t) / (max_val - min_val)

    for idx, val in enumerate(arr):
        arr[idx] = (val - min_val) * ratio + min_t
    return arr


def init_directed_network(file_dir, top_n=10):
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

    f_n_path_name = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_path_prob = os.path.join(file_dir, "output", "pathway_prob.csv")

    p_n = np.genfromtxt(f_n_path_name, dtype=str,
                        delimiter=',', max_rows=top_n)
    p_p = np.genfromtxt(f_n_path_prob, dtype=float,
                        delimiter=',', max_rows=top_n)

    d_f = pd.DataFrame(np.transpose([p_n, p_p]), columns=['name', 'prob'])

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

        matched_spe = re.findall(r"S(\d+)", val['name'])
        matched_reaction = re.findall(r"R(\d+)", val['name'])

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


if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)
    TOP_N = 100

    RN_OBJ = init_directed_network(FILE_DIR, top_n=TOP_N)
    network_to_gephi_input_file(
        RN_OBJ, FILE_DIR, "network_" + str(TOP_N) + ".gexf")

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
