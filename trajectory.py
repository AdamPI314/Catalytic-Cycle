"""
get trajectory related information
"""

import os
import sys
from collections import defaultdict, OrderedDict
import numpy as np
from copy import deepcopy
from scipy.interpolate import CubicSpline
import scipy.optimize as opt
import parse_spe_reaction_info as psri
import interpolation
import global_settings


def convert_path_prob_to_concentration(data_dir, atom_followed="C",
                                       path_prob=None, default_coef=None):
    """
    if default_coef is not None, use it as default coefficient
    convert total pathway probability to concentration
    for example, C3H8, suppose [C3H8] = 1.0 and we are following "C"
    atom, then the corresponding total pathway probability should be
    1.0 * 3, since each C3H8 has 3 "C" atoms, in other word, concentration
    should be pathway probability divide by 3.0
    """
    if path_prob is None:
        return None
    if path_prob is []:
        return None

    _, spe_n_i_d = psri.parse_spe_info(data_dir)
    spe_composition = psri.read_spe_composition(
        os.path.join(data_dir, "input", "spe_composition.json"))

    spe_idx_coefficient = dict()
    for _, val in enumerate(spe_composition):
        if atom_followed in spe_composition[val]:
            spe_idx_coefficient[spe_n_i_d[val]] = float(
                spe_composition[val][atom_followed])
        else:
            spe_idx_coefficient[spe_n_i_d[val]] = 0.0

    if default_coef is not None:
        for val in spe_idx_coefficient:
            if spe_idx_coefficient[val] != 0:
                spe_idx_coefficient[val] = default_coef

    if np.shape(path_prob)[0] > 0:
        if np.shape(path_prob[0]) is ():
            print("1D array", "shape:\t", len(path_prob))
            for idx, _ in enumerate(path_prob):
                if float(spe_idx_coefficient[str(idx)]) != 0:
                    path_prob[idx] /= float(spe_idx_coefficient[str(idx)])

    return path_prob


def convert_concentration_to_path_prob(data_dir, atom_followed="C", spe_conc=None,
                                       renormalization=True, default_coef=None):
    """
    if default_coef is not None, use it as default coefficient
    convert concentration to corresponding total pathway probability
    for example, C3H8, suppose [C3H8] = 1.0 and we are following "C"
    atom, then the corresponding total pathway probability should be
    1.0 * 3, since each C3H8 has 3 "C" atoms
    Warning: spe_conc should be read from dlsode calculation, it is
    guaranteed outside that dimensions of spe_conc match the mechanism
    """
    if spe_conc is None:
        return None
    if spe_conc is []:
        return None

    path_prob = deepcopy(spe_conc)

    _, spe_n_i_d = psri.parse_spe_info(data_dir)
    spe_composition = psri.read_spe_composition(
        os.path.join(data_dir, "input", "spe_composition.json"))

    spe_idx_coefficient = dict()
    for _, val in enumerate(spe_composition):
        if atom_followed in spe_composition[val]:
            spe_idx_coefficient[spe_n_i_d[val]] = float(
                spe_composition[val][atom_followed])
        else:
            spe_idx_coefficient[spe_n_i_d[val]] = 0.0
    # print(spe_composition, spe_idx_coefficient)

    if default_coef is not None:
        for val in spe_idx_coefficient:
            if spe_idx_coefficient[val] != 0:
                spe_idx_coefficient[val] = default_coef
    # print(spe_composition, spe_idx_coefficient)

    if np.shape(path_prob)[0] > 0:
        if np.shape(path_prob[0]) is ():
            print("1D array", "shape:\t", len(path_prob))
            for idx, _ in enumerate(path_prob):
                if float(spe_idx_coefficient[str(idx)]) != 0:
                    path_prob[idx] *= float(spe_idx_coefficient[str(idx)])
                else:
                    path_prob[idx] *= 0.0
            if renormalization is True:
                path_prob /= np.sum(path_prob)
        else:
            print("2D array", "shape:\t", np.shape(path_prob))
            for idx in range(np.shape(path_prob)[1]):
                if float(spe_idx_coefficient[str(idx)]) != 0:
                    path_prob[:, idx] *= float(spe_idx_coefficient[str(idx)])
                else:
                    path_prob[:, idx] *= 0.0
            if renormalization is True:
                for idx, _ in enumerate(path_prob):
                    path_prob[idx, :] /= np.sum(path_prob[idx, :])

    return path_prob


def get_species_with_top_n_concentration(data_dir, exclude, top_n=10, traj_max_t=100.0,
                                         tau=10.0, end_t=1.0, tag="M", atoms=None):
    """
    get species concentration at a tau, where tau is the ratio of the time_wanted/end_time
    """
    if atoms is None:
        atoms = ["C"]
    if exclude is None:
        exclude = []

    time = np.loadtxt(os.path.join(data_dir, "output",
                                   "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    conc_all = np.loadtxt(os.path.join(data_dir, "output",
                                       "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")

    n_spe = np.shape(conc_all)[1]

    data = [float] * n_spe
    for i in range(n_spe):
        data[i] = interpolation.interp1d(time, conc_all[:, i], tau * end_t)

    c_idx_map = defaultdict(set)
    for idx, val in enumerate(data):
        c_idx_map[val].add(str(idx))
    c_idx_map = OrderedDict(sorted(c_idx_map.items(), reverse=True))

    spe_idx_name_dict, _ = psri.parse_spe_info(data_dir)
    spe_composition = psri.read_spe_composition(
        os.path.join(data_dir, "input", "spe_composition.json"))

    spe_idx_list = []
    counter = 0

    for _, val in enumerate(c_idx_map):
        if counter < top_n:
            spe_idx = next(iter(c_idx_map[val]))
            indicator = False
            for _, atom in enumerate(atoms):
                if atom in spe_composition[spe_idx_name_dict[spe_idx]]:
                    indicator = True
                    break
            if spe_idx_name_dict[spe_idx] not in exclude and indicator:
                print(val, spe_idx, spe_idx_name_dict[spe_idx])
                spe_idx_list.append(int(spe_idx))
                counter += 1

    # species doesn't contain atom we are interested in
    exclude_spe_name_list = []
    for idx, s_n_t in enumerate(spe_composition):
        indicator = False
        for _, atom in enumerate(atoms):
            if atom in spe_composition[s_n_t]:
                indicator = True
        if indicator is False:
            exclude_spe_name_list.append(s_n_t)
    spe_name_list = [str(spe_idx_name_dict[str(x)]) for x in spe_idx_list]
    return spe_idx_list, spe_name_list, exclude_spe_name_list


def get_normalized_concentration(data_dir, tag="fraction", exclude_names=None, renormalization=True):
    """
    return normalized concentration
    """
    if exclude_names is None:
        exclude_names = []
    conc = np.loadtxt(os.path.join(data_dir, "output",
                                   "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")
    if renormalization is False:
        return conc

    _, s_n_idx = psri.parse_spe_info(data_dir)
    # renormalization
    exclude_idx_list = [int(s_n_idx[x]) for x in exclude_names]
    # set the concentration of these species to be zero
    for _, idx in enumerate(exclude_idx_list):
        conc[:, idx] = 0.0
    # normalize
    for idx, _ in enumerate(conc):
        conc[idx, :] /= np.sum(conc[idx, :])

    return conc


def get_normalized_concentration_at_time(data_dir, tag="fraction", tau=10.0, end_t=1.0, exclude_names=None, renormalization=True):
    """
    return normalized species concentration at time
    """
    if exclude_names is None:
        exclude_names = []
    time = np.loadtxt(os.path.join(data_dir, "output",
                                   "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    conc_all = np.loadtxt(os.path.join(data_dir, "output",
                                       "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")

    n_spe = np.shape(conc_all)[1]

    conc = [float] * n_spe
    for i in range(n_spe):
        conc[i] = interpolation.interp1d(time, conc_all[:, i], tau * end_t)

    _, s_n_idx = psri.parse_spe_info(data_dir)
    exclude_idx_list = [int(s_n_idx[x]) for x in exclude_names]
    # set the concentration of these species to be zero
    for _, idx in enumerate(exclude_idx_list):
        conc[idx] = 0.0

    if renormalization is False:
        return conc

    _, s_n_idx = psri.parse_spe_info(data_dir)
    # renormalization
    exclude_idx_list = [int(s_n_idx[x]) for x in exclude_names]
    # set the concentration of these species to be zero
    for _, idx in enumerate(exclude_idx_list):
        conc[idx] = 0.0
    conc /= np.sum(conc)
    return conc


def get_time_at_temperature_differential_maximum(data_dir, l_b=0.7, h_b=0.8):
    """
    return time at which the first order differential of temperature is maximum
    l_b: lower bound
    h_b: higher bound
    """
    print(data_dir)
    f_n_time = os.path.join(data_dir, "output", "time_dlsode_M.csv")
    f_n_temp = os.path.join(data_dir, "output", "temperature_dlsode_M.csv")

    time = np.loadtxt(f_n_time, dtype=float, delimiter=',')
    temp = np.loadtxt(f_n_temp, dtype=float, delimiter=',')
    print(np.shape(time))

    temp_gradient = np.gradient(temp, time)
    print(temp_gradient)

    # cubic spline function
    func = CubicSpline(time, temp_gradient)

    # maximum, multiply by -1, maximum --> minimum
    # multivariate method
    # https://docs.scipy.org/doc/scipy/reference/optimize.html
    # be careful when setting the bounds, avoid stiffness problem
    max_x = opt.fmin_l_bfgs_b(lambda x: -1 * func(x), time[-1],
                              bounds=[(l_b, h_b)], approx_grad=True)
    print(max_x)

    if "CONVERGENCE" in str(max_x[-1]['task']):
        print("convergence achieved\t", max_x[0][0])
        return max_x[0][0]
    else:
        print("not converges\t")


if __name__ == '__main__':
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    G_S = global_settings.get_setting(DATA_DIR)
    print(DATA_DIR)
    # get_normalized_concentration_at_time(
    #     DATA_DIR, tag="M", end_t=0.9, exclude_names=None, renormalization=True)
    # convert_concentration_to_path_prob(
    #     DATA_DIR, atom_followed="C", spe_conc=[1.0, 2.0], renormalization=True)
    get_time_at_temperature_differential_maximum(DATA_DIR)
