"""
parse species patten, given a path
"""

import os
import sys
import re
from itertools import combinations


def parse_path_length(path):
    """
    parse path length
    """
    matched_tmp = re.findall(r"(S\d+)", path)
    return len(matched_tmp)


def parse_species(path):
    """
    parse sepcies name, return a dictionary of species and count
    """
    matched_tmp = re.findall(r"(S\d+)", path)
    d_map = dict()
    # not the first and last element
    for _, val in enumerate(matched_tmp[1:-1]):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_species_production_reaction(path, spe):
    """
    parse sepcies name production, return a dictionary of species and count
    reactions ends with spe
    """
    matched_tmp = re.findall(r"(R\d+)" + spe + r"(?=$|R)", path)
    d_map = dict()
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_species_production_path(path, spe):
    """
    parse sepcies name production, return a dictionary of species and count
    pathway or sub-pathway ends with spe
    """
    idx_tmp = [(m.start(0), m.end(0))
               for m in re.finditer(spe + r"(?=$|R)", path)]
    d_map = dict()
    # from intitial species to that species
    for _, val in enumerate(idx_tmp):
        sub_path = path[0:val[1]]
        if sub_path not in d_map:
            d_map[sub_path] = 1
        else:
            d_map[sub_path] += 1
    return d_map


def parse_reaction(path):
    """
    parse reactions, return a dictionary of reactions and their count
    """
    matched_tmp = re.findall(r"(R\d+)", path)
    d_map = dict()
    # include all reactions
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    return d_map


def parse_initiation_reaction(path):
    """
    parse initiation_reactions, return a dictionary of reactions and their count
    """
    matched_tmp = re.findall(r"(R\d+)", path)
    d_map = dict()
    # include all reactions
    if len(matched_tmp) >= 1:
        for _, val in enumerate(matched_tmp[0:1]):
            if val not in d_map:
                d_map[val] = 1
            else:
                d_map[val] += 1
    return d_map


def parse_species_cycle(path):
    """
    parse sepcies name, return a dictionary of species cycle and count
    """
    matched_tmp = re.findall(r"(S\d+)", path)
    d_map = dict()
    for _, val in enumerate(matched_tmp):
        if val not in d_map:
            d_map[val] = 1
        else:
            d_map[val] += 1
    cycle_map = dict()
    for key, val in d_map.items():
        # search pathway if occurence of spcies >= 2
        if val >= 2:
            idx_tmp = [(m.start(0), m.end(0)) for m in re.finditer(key, path)]
            for p_x in combinations(idx_tmp, 2):
                cycle = path[p_x[0][0]:p_x[1][1]]
                if cycle not in cycle_map:
                    cycle_map[cycle] = 1
                else:
                    cycle_map[cycle] += 1

    return cycle_map


if __name__ == "__main__":
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(parse_path_length("S10"))
