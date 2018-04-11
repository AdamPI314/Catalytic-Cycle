"""
parse species patten, given a path
"""

import os
import sys
import re
from itertools import combinations
import parse_spe_reaction_info as psri
import atom_scheme as asch


def path_contain_regex(path, path_reg=None):
    """
    test if a path string contains regular expression 
    """
    if path_reg is None:
        return True
    return bool(re.search(path_reg, path))


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


def get_spe_production_sub_path(pathname="S60R-100001S90R1162S94", net_r=None, net_p=None, spe_idx=10, s_p_r_c=None, n_threshold=0):
    """
    return sub-pathway that give rise to a species, produces a species along a pathway
    """
    if net_r is None or net_p is None:
        return []

    # get rid of R-1000003S90, don't need it here
    pathname = re.sub(r"R-\d+S\d+", r'', pathname)

    sub_path = []
    sub_path_reaction = []
    # parse pathway
    p = re.compile(r"R(\d+)S(\d+)")

    reaction_str = ""
    for m in p.finditer(pathname):
        # print(m.start(), m.end(), m.group(), m.group(1))
        r_idx = m.group(1)
        reaction_str += ("R" + str(r_idx))
        path_tmp = str(pathname[0:m.end()])

        # gotta to make sure the last reaction creates at least one OH
        if str(spe_idx) in net_p[str(r_idx)]:
            # the whole sub-pathway create net number of OH more than threshold
            if parse_net_species_along_path_using_reaction(
                    path_tmp, net_r=net_r, net_p=net_p, spe_idx=spe_idx, s_p_r_c=s_p_r_c) >= n_threshold:
                sub_path.append(path_tmp)
                sub_path_reaction.append(reaction_str)

    # print(sub_path, sub_path_reaction)
    return sub_path, sub_path_reaction


def parse_species_along_path_using_reaction(pathname="S60R-100001S90R1162S94", net_r_p=None, spe_idx=10, s_p_r_c=None):
    """
    calculate number of species being used or produced along a path, depends on net_r_p,
    if net_r_p is net_reactant, return net species be consumed
    if net_r_p is net_product, return net species be produced
    notice OH might be not directly shown on a path, but can be side products of reactions from path
    """
    if net_r_p is None:
        return 0

    number = 0

    # chattering might produce desired spe_idx as well
    # match S1R-1000S2 combination
    m_srs = re.findall(r"(S\d+R-\d+S\d+)", pathname)
    for _, s_r_s in enumerate(m_srs):
        # print(s_r_s)
        s_1 = next(iter(re.findall(r"S(\d+)R-\d+S\d+", s_r_s)), 0)
        s_2 = next(iter(re.findall(r"S\d+R-\d+S(\d+)", s_r_s)), 0)

        if s_1 != s_2 and (s_1, s_2) in s_p_r_c:
            for pair_idx in s_p_r_c[(s_1, s_2)]:
                r_idx_c = s_p_r_c[(s_1, s_2)][pair_idx]['r_idx']
                # print(r_idx_c)
                if str(spe_idx) in net_r_p[r_idx_c]:
                    # print("bingo")
                    number += int(net_r_p[r_idx_c][str(spe_idx)])

    # get rid of R-1000003S90, don't need it here
    pathname = re.sub(r"R-\d+S\d+", r'', pathname)

    # parse pathway
    matched_reaction = re.findall(r"R(\d+)", pathname)

    for _, r_idx in enumerate(matched_reaction):
        print(r_idx)
        if str(r_idx) in net_r_p:
            # print(net_r_p[r_idx])
            if str(spe_idx) in net_r_p[str(r_idx)]:
                # print([net_r_p[r_idx][str(spe_idx)]])
                number += int(net_r_p[r_idx][str(spe_idx)])

    return number


def parse_net_species_along_path_using_reaction(pathname="S60R-100001S90R1162S94", net_r=None, net_p=None, spe_idx=10, s_p_r_c=None):
    """
    calculate net number of species being produced along a path
    net_r is net_reactant
    net_p is net_product
    notice OH might be not directly shown on a path, but can be side products of reactions from path
    """
    if net_r is None or net_p is None:
        return 0

    number = 0

    # chattering might produce desired spe_idx as well
    # match S1R-1000S2 combination
    m_srs = re.findall(r"(S\d+R-\d+S\d+)", pathname)
    for _, s_r_s in enumerate(m_srs):
        # print(s_r_s)
        s_1 = next(iter(re.findall(r"S(\d+)R-\d+S\d+", s_r_s)), 0)
        s_2 = next(iter(re.findall(r"S\d+R-\d+S(\d+)", s_r_s)), 0)

        if s_1 != s_2 and (s_1, s_2) in s_p_r_c:
            for pair_idx in s_p_r_c[(s_1, s_2)]:
                r_idx_c = s_p_r_c[(s_1, s_2)][pair_idx]['r_idx']
                # print(r_idx_c)
                if str(spe_idx) in net_r[r_idx_c]:
                    number -= int(net_r[r_idx_c][str(spe_idx)])
                if str(spe_idx) in net_p[r_idx_c]:
                    number += int(net_p[r_idx_c][str(spe_idx)])

    # get rid of R-1000003S90, don't need it here
    pathname = re.sub(r"R-\d+S\d+", r'', pathname)

    # parse pathway
    matched_reaction = re.findall(r"R(\d+)", pathname)

    for _, r_idx in enumerate(matched_reaction):
        # print(r_idx)
        if str(r_idx) in net_r:
            if str(spe_idx) in net_r[str(r_idx)]:
                number -= int(net_r[r_idx][str(spe_idx)])
        if str(r_idx) in net_p:
            if str(spe_idx) in net_p[str(r_idx)]:
                number += int(net_p[r_idx][str(spe_idx)])

    return number


def calculate_reaction_branching_number(r_idx=0, net_reactant=None, net_product=None, s_idx_name=None, atom_scheme=None, atom_followed="C"):
    """
    Assume reaction occurs, thereafter there is nothing to do with reaction rate
    R->S1 + S2, the number for RS1 will be 1+1 = 2, assuming S1 and S2 both contain only one atom being followed
    R->2S1 + S2, the number for RS1 will be  2+1=3, assuming S1 and S2 both contain only one atom being followed
    A important thing is a net product can have at most one followed atom otherwise it it meanless to do it
    since the idea here is branching, 1->1, or 1->2 transition
    """
    if net_reactant is None:
        return 0
    if net_product is None:
        return 0
    if s_idx_name is None:
        return 0
    if atom_scheme is None:
        return 0
    if atom_followed not in atom_scheme:
        return 0

    n_r = net_reactant[str(r_idx)]
    n_p = net_product[str(r_idx)]

    num_r = 0
    num_p = 0
    for _, s_idx in enumerate(n_r):
        s_name = s_idx_name[str(s_idx)]
        if str(s_name) in atom_scheme[atom_followed]:
            if int(atom_scheme[atom_followed][str(s_name)]) >= 1:
                print(s_name, atom_scheme[atom_followed][str(s_name)])
                # num_r += 1
                num_r += float(n_r[s_idx])
    for _, s_idx in enumerate(n_p):
        s_name = s_idx_name[str(s_idx)]
        if str(s_name) in atom_scheme[atom_followed]:
            if int(atom_scheme[atom_followed][str(s_name)]) >= 1:
                print(s_name, atom_scheme[atom_followed][str(s_name)])
                # num_p += 1
                num_p += float(n_p[s_idx])

    if num_r == 0:
        return 0
    return float(num_p) / float(num_r)


def calculate_path_branching_number(pathname="S60R-100001S90R1162S94", net_reactant=None,
                                    net_product=None, s_idx_name=None, atom_scheme=None, atom_followed="C"):
    """
    calculate path branching number
    """
    if net_reactant is None:
        return 1.0
    if net_product is None:
        return 1.0
    if s_idx_name is None:
        return 1.0
    if atom_scheme is None:
        return 1.0
    if atom_followed not in atom_scheme:
        return 1.0

    number = 1.0
    # get rid of R-1000003S90, don't need it here
    pathname = re.sub(r"R-\d+S\d+", r'', pathname)

    # parse pathway
    matched_reaction = re.findall(r"R(\d+)", pathname)

    for _, r_idx in enumerate(matched_reaction):
        print(r_idx)
        number *= calculate_reaction_branching_number(r_idx=r_idx, net_reactant=net_reactant, net_product=net_product,
                                                      s_idx_name=s_idx_name, atom_scheme=atom_scheme, atom_followed=atom_followed)

    return number


if __name__ == "__main__":
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    # print(parse_path_length("S10"))
    # NET_REACTANT = psri.parse_reaction_net_reactant(DATA_DIR)
    # NET_PRODUCT = psri.parse_reaction_net_product(DATA_DIR)
    # ATOM_SCHEME = asch.get_atom_scheme(DATA_DIR)
    # S_IDX_NAME, _ = psri.parse_spe_info(DATA_DIR)
    # S_P_R_C = psri.parse_species_pair_reaction(DATA_DIR)
    # parse_spe_production_along_path(net_product=NET_PRODUCT)

    # calculate_reaction_branching_number(
    #     r_idx=452, net_reactant=NET_REACTANT, net_product=NET_PRODUCT,
    #     s_idx_name=S_IDX_NAME, atom_scheme=ATOM_SCHEME, atom_followed="C")

    # calculate_path_branching_number(
    #     pathname="S60R-100001S90R1162S94R-100006S101R1222S46R452S17", net_reactant=NET_REACTANT, net_product=NET_PRODUCT,
    #     s_idx_name=S_IDX_NAME, atom_scheme=ATOM_SCHEME, atom_followed="C")

    # parse_species_along_path_using_reaction(
    #     pathname="S90R1162S94R-100006S101R1222S46R90S14", net_r_p=NET_REACTANT, spe_idx=10, s_p_r_c=S_P_R_C)

    # get_spe_production_sub_path(
    #     pathname="S90R1162S94R-100006S101R1222S46R90S14R90S17R90S45", net_p=NET_REACTANT, spe_idx=10, s_p_r_c=S_P_R_C)

    print(path_contain_regex('S60R1162S94', 'S(39|50)'))
