"""
network related, generate network, save to json, save to gephi compatible file
"""

import os
import time
import sys
import re
import numpy as np
import pandas as pd
import networkx as nx


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


def init_network(file_dir, top_n=10):
    """
    init network
    return networkx.MultiDiGraph
    """
    f_n_path_name = os.path.join(file_dir, "output", "pathway_name.csv")
    f_n_path_prob = os.path.join(file_dir, "output", "pathway_prob.csv")

    p_n = np.genfromtxt(f_n_path_name, dtype=str,
                        delimiter=',', max_rows=top_n)
    p_p = np.genfromtxt(f_n_path_prob, dtype=float,
                        delimiter=',', max_rows=top_n)

    d_f = pd.DataFrame(np.transpose([p_n, p_p]), columns=['name', 'prob'])
    # print(d_f)

    # directed graph
    d_g = nx.MultiDiGraph()

    # record all nodes
    nodes = set()
    for _, val in d_f.iterrows():
        matched_spe = re.findall(r"S(\d+)", val['name'])
        for _, spe in enumerate(matched_spe):
            nodes.add(change_spe_name(spe))

    for _, val in enumerate(nodes):
        d_g.add_node(val, weight=0.0, label=str(val))

    for _, val in d_f.iterrows():
        prob = float(val['prob'])

        matched_spe = re.findall(r"S(\d+)", val['name'])
        matched_reaction = re.findall(r"R(\d+)", val['name'])

        for idx, spe in enumerate(matched_spe):
            d_g.node[change_spe_name(spe)]['weight'] += 1.0 * prob
            if idx > 0:
                src = change_spe_name(matched_spe[idx - 1])
                des = change_spe_name(spe)
                rxn = change_rxn_name(matched_reaction[idx - 1])
                if d_g.has_edge(src, des):
                    if rxn in d_g[src][des].keys():
                        # print("found")
                        d_g[src][des][rxn]['weight'] += 1.0 * prob
                    else:
                        d_g[src][des].update(
                            {rxn: {'Label': rxn, 'weight': 1.0 * prob}})
                else:
                    d_g.add_edge(src, des, key=rxn,
                                 Label=rxn, weight=1.0 * prob)

    # print(d_g.edges())
    # for idx, val in enumerate(d_g.nodes()):
    #     print(idx, val, d_g.node[val]['weight'])
    # for idx, val in enumerate(d_g.edges()):
    #     print(idx, val, d_g[val[0]][val[1]].items())

    return d_g


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

    nx_obj = init_network(FILE_DIR, top_n=10)
    network_to_gephi_input_file(nx_obj, FILE_DIR, "network.gexf")

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
