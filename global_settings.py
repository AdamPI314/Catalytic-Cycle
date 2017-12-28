"""
global settings
"""
import sys
import os
from collections import OrderedDict, defaultdict
import union_find


def get_fast_rxn_trapped_spe(file_dir):
    """
    get_fast_rxn_trapped_spe
    """
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import local_settings
        return local_settings.get_fast_rxn_trapped_spe()
    except IOError:
        return OrderedDict(), OrderedDict()


def get_union_find_group(file_dir):
    """
    return union_find_groups
    """
    _, trapped_species = get_fast_rxn_trapped_spe(file_dir)

    counter = 0
    spe_idx_label = dict()
    label_spe_idx = dict()

    u_set = set()

    for _, val in enumerate(trapped_species):
        if int(val) not in u_set:
            u_set.add(int(val))
            spe_idx_label[int(val)] = counter
            label_spe_idx[counter] = int(val)
            counter += 1
        if trapped_species[val] not in u_set:
            u_set.add(int(trapped_species[val]))
            spe_idx_label[int(trapped_species[val])] = counter
            label_spe_idx[counter] = int(trapped_species[val])
            counter += 1
    # print(spe_idx_label, label_spe_idx)
    wqnpc = union_find.WeightedQuickUnionWithPathCompression(len(u_set))
    for _, val in enumerate(trapped_species):
        idx1 = int(spe_idx_label[int(val)])
        idx2 = int(spe_idx_label[int(trapped_species[val])])
        wqnpc.unite(idx1, idx2)

    # unique labels
    unique_labels = set()
    for idx, _ in enumerate(spe_idx_label):
        l_tmp = wqnpc.root(idx)
        unique_labels.add(l_tmp)

    # unique labels and their group
    unique_labels_group = defaultdict(set)
    for idx, _ in enumerate(spe_idx_label):
        l_tmp = wqnpc.root(idx)
        if l_tmp not in unique_labels_group:
            unique_labels_group[l_tmp] = set()
            unique_labels_group[l_tmp].add(label_spe_idx[idx])
        else:
            unique_labels_group[l_tmp].add(label_spe_idx[idx])
    # print(unique_labels_group)

    # species index and the big group it belongs to
    idx_group = defaultdict(set)
    for idx, _ in enumerate(spe_idx_label):
        if idx in unique_labels_group:
            idx_group[str(label_spe_idx[idx])] = unique_labels_group[idx]
        else:
            idx_group[str(label_spe_idx[idx])
                      ] = unique_labels_group[int(wqnpc.root(int(idx)))]
    print(idx_group)

    return idx_group


def get_setting(file_dir):
    """
    return global settings
    """
    setting = {}
    try:
        sys.path.append(os.path.join(file_dir, "input"))
        import local_settings
        return local_settings.get_local_settings()
    except Exception:
        return setting
