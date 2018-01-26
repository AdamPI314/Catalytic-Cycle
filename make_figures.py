"""
plot jobs
"""
import os
import sys
import re
from copy import deepcopy
import numpy as np
import json
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt

import parse_spe_reaction_info as psri
import trajectory
import pattern_statistics
import global_settings
from tools import get_colors_markers_linestyles, pathway_time_2_array_index
import naming
import interpolation


def plot_path_length_statistics(file_dir, init_spe=62, atom_followed="C",
                                end_t=1.0, end_spe=None):
    """
    plot path length statistics
    """
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    if suffix is not None:
        suffix += "_S" + str(end_spe)
    in_f_n = os.path.join(file_dir, "output", "path_length" + suffix + ".csv")
    fig_name = os.path.join(
        file_dir, "output", "path_length" + suffix + ".jpg")

    colors, markers, _ = get_colors_markers_linestyles()

    data = np.loadtxt(in_f_n, dtype=float, delimiter=',')

    data_x = data[:, 0]
    data_y = data[:, 1]
    # data_y /= np.sum(data_y)

    spe_idx_n_dict, _ = psri.parse_spe_info(file_dir)
    # specify label for lines
    labels = [spe_idx_n_dict[str(end_spe)]]

    delta_n = int(len(data_x) / 50)
    if delta_n is 0:
        delta_n = 1

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    idx = 0
    a_x.plot(data_x[::delta_n], data_y[::delta_n],
             color=colors[idx], marker=markers[idx], label=labels[idx])

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 10.0})
    leg.get_frame().set_alpha(0.7)

    a_x.set_xlim([data_x[0], data_x[-1]])
    a_x.grid()

    a_x.set_xlabel("path length measured by number of species")
    a_x.set_ylabel("frequency")
    a_x.set_title("path length distribution")

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, fig_name), dpi=500)
    plt.close()

    return


def plot_cumulative_pathway_prob(file_dir, init_spe=62, atom_followed="C", end_t=1.0, tau=10.0, species_path=True, top_n=10, time_axis=0):
    """
    plot cumulative pathway probability at a time point
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)

    f_n_pp = os.path.join(file_dir, "output",
                          prefix + "pathway_prob" + suffix + ".csv")
    data = np.genfromtxt(f_n_pp, delimiter=",", dtype=float)

    # data dimensions, could be array or matrix
    dim_n = len(np.shape(data))
    if dim_n == 1:
        data = data
    elif dim_n == 2:
        data = data[:, time_axis]
    else:
        return

    data = sorted(data, reverse=True)[0:top_n]

    data_c = deepcopy(data)
    data_c[0] = 0.0
    for idx, _ in enumerate(data_c):
        if idx >= 1:
            data_c[idx] += data_c[idx - 1]
    print(data_c[-1])

    fig, a_x = plt.subplots(1, 1, sharex=False, sharey=False)

    end_point = int(end_t * len(data))
    a_x.plot(data[1:end_point:1], "-.")
    a_x.plot(data_c[1:end_point:1], "-*")
    a_x.grid()

    suffix += "_top_" + str(top_n)
    f_n_out = os.path.join(file_dir, "output", prefix +
                           "pathway_probability" + suffix + ".jpg")
    fig.savefig(f_n_out, dpi=500)

    return


def plot_concentrations(file_dir, spe_idx=None, tau=10.0, end_t=1.0, tag="fraction", exclude_names=None,
                        renormalization=True, semilogy=False, hasTemp=True):
    """
    plot concentrations give species index list, if exclude is not None, means we are going
    to renormalize the molelar fraction
    """
    if exclude_names is None:
        exclude_names = []

    spe_idx_tmp = deepcopy(spe_idx)
    if spe_idx_tmp is None:
        spe_idx_tmp = [0]

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)

    if hasTemp is True:
        s_idx_n["-1"] = "Temp"
        spe_idx_tmp.append(-1)

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    conc = trajectory.get_normalized_concentration(
        file_dir, tag=tag, exclude_names=exclude_names, renormalization=renormalization)

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx_tmp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            if semilogy is True:
                a_x_left.semilogy(time[0:end_point], conc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            else:
                a_x_left.plot(time[0:end_point], conc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                              color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.set_xlim([0, tau * end_t])

    a_x_left.set_xlabel("Time/sec")
    a_x_left.set_ylabel("[X]")

    if hasTemp is True:
        leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
        leg_right.get_frame().set_alpha(0.7)
        a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx_tmp)
    # plt.title(s_n_str)

    fig.savefig(os.path.join(file_dir, "output",
                             "trajectory_" + s_n_str + ".jpg"), dpi=500)
    plt.close()


def plot_spe_concentrations_derivative(file_dir, spe_idx=None, tau=10.0, end_t=1.0,
                                       tag="fraction", exclude_names=None, renormalization=True):
    """
    plot concentrations give species index list, if exclude is not None, means we are going
    to renormalize the molelar fraction
    """
    if exclude_names is None:
        exclude_names = []

    spe_idx_tmp = deepcopy(spe_idx)
    if spe_idx_tmp is None:
        spe_idx_tmp = [0]

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)
    s_idx_n["-1"] = "Temp"

    spe_idx_tmp.append(-1)

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    conc = trajectory.get_normalized_concentration(
        file_dir, tag=tag, exclude_names=exclude_names, renormalization=renormalization)

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx_tmp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]

            data_y = np.gradient(conc[:, s_idx], time)
            # a_x_left.semilogy(time[0:end_point], data_y[0:end_point], marker=m_k, markevery=delta_n,
            #                   color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            a_x_left.plot(time[0:end_point], data_y[0:end_point], marker=m_k, markevery=delta_n,
                          color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    leg_left = a_x_left.legend(loc=9, fancybox=True, prop={'size': 10.0})
    leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()

    a_x_left.set_xlabel("Time/sec")

    a_x_left.set_ylabel("[$\dot X$]")
    a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx_tmp)
    # plt.title(s_n_str)

    fig.savefig(os.path.join(file_dir, "output",
                             "spe_concentrations_derivative_" + s_n_str + ".jpg"), dpi=500)
    plt.close()


def plot_spe_drc(file_dir, spe_idx=None, tau=10.0, end_t=1.0, tag="fraction", reciprocal=False):
    """
    plot species destruction rate constant, give species index list
    """
    spe_idx_tmp = deepcopy(spe_idx)
    if spe_idx_tmp is None:
        spe_idx_tmp = [0]

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)
    s_idx_n["-1"] = "Temp"

    spe_idx_tmp.append(-1)

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    spe_drc = np.loadtxt(os.path.join(file_dir, "output",
                                      "drc_dlsode_" + str(tag) + ".csv"), delimiter=",")
    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx_tmp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point], markevery=delta_n,
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            if reciprocal is False:
                a_x_left.semilogy(time[0:end_point], spe_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            else:
                a_x_left.semilogy(time[0:end_point], 1.0 / spe_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    if reciprocal is False:
        leg_left = a_x_left.legend(loc=9, fancybox=True, prop={'size': 10.0})
    else:
        leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})

    leg_right = a_x_right.legend(loc=4, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.set_xlim([0.05, time[end_point]])

    a_x_left.set_xlabel("time/s")
    if reciprocal is False:
        a_x_left.set_ylabel("k/s$^{-1}$")
    else:
        a_x_left.set_ylabel("k$^{-1}/s$")

    a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx_tmp)
    # plt.title(s_n_str)

    if reciprocal is False:
        fig.savefig(os.path.join(file_dir, "output",
                                 "spe_drc_" + s_n_str + ".jpg"), dpi=500)
    else:
        fig.savefig(os.path.join(file_dir, "output",
                                 "spe_drc_reciprocal_" + s_n_str + ".jpg"), dpi=500)

    plt.close()


def plot_chattering_group_drc(file_dir, tau=10.0, end_t=1.0, tag="fraction", reciprocal=False):
    """
    plot chattering group destruction rate constant
    """
    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(file_dir)
    with open(os.path.join(file_dir, "output", "chattering_group_info.json"), 'r') as f_h:
        chattering_group_info = json.load(f_h)

    chattering_group_idx_n = dict()
    for idx1, val1 in enumerate(chattering_group_info):
        label_tmp = ""
        for idx2, val2 in enumerate(chattering_group_info[val1]):
            s_idx_tmp = chattering_group_info[val1][val2]
            label_tmp += s_idx_n[s_idx_tmp]
            if idx2 < len(chattering_group_info[val1]) - 1:
                label_tmp += "+"
        chattering_group_idx_n[str(idx1)] = label_tmp
    chattering_group_idx_n["-1"] = "Temp"

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    c_g_drc = np.loadtxt(os.path.join(file_dir, "output",
                                      "chattering_group_drc_dlsode_" + str(tag) + ".csv"), delimiter=",")

    group_id_temp = [i for i in range(np.shape(c_g_drc)[1])]
    group_id_temp.append(-1)

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 10)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in group_id_temp:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=chattering_group_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            if reciprocal is False:
                a_x_left.semilogy(time[0:end_point], c_g_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=chattering_group_idx_n[str(s_idx)])
            else:
                a_x_left.semilogy(time[0:end_point], 1.0 / c_g_drc[0:end_point, s_idx], marker=m_k, markevery=delta_n,
                                  color=colors[counter % (len(colors) - 1)], label=chattering_group_idx_n[str(s_idx)])

            counter += 1
    if reciprocal is False:
        leg_left = a_x_left.legend(loc=9, fancybox=True, prop={'size': 10.0})
    else:
        leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})

    leg_right = a_x_right.legend(loc=4, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()
    a_x_left.set_xlim([0.05, time[end_point]])

    a_x_left.set_xlabel("time/s")

    if reciprocal is False:
        a_x_left.set_ylabel("k/s$^{-1}$")
    else:
        a_x_left.set_ylabel("k$^{-1}$/s")

    a_x_right.set_ylabel("T/K")

    if reciprocal is False:
        fig.savefig(os.path.join(file_dir, "output",
                                 "chattering_group_drc" + ".jpg"), dpi=500)
    else:
        fig.savefig(os.path.join(file_dir, "output",
                                 "chattering_group_drc_reciprocal" + ".jpg"), dpi=500)

    plt.close()


def plot_reaction_rates(file_dir, reaction_idx=None, tau=10.0, end_t=1.0, tag="fraction"):
    """
    plot reaction rates give reaction index list
    """

    colors, markers, _ = get_colors_markers_linestyles()

    _, rxn_idx_n = psri.parse_reaction_and_its_index(file_dir)
    rxn_idx_n["-1"] = "Temp"
    reaction_idx.append(-1)

    if reaction_idx is None:
        reaction_idx = [0]
    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    rxn_rates = np.loadtxt(os.path.join(file_dir, "output",
                                        "reaction_rate_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    counter = 0
    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    delta_n = int(end_point / 25)
    if delta_n is 0:
        delta_n = 1

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in reaction_idx:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point:delta_n], temp[0:end_point:delta_n],
                           color=colors[-1], label=rxn_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            a_x_left.semilogy(time[0:end_point:delta_n], rxn_rates[0:end_point:delta_n, s_idx], marker=m_k,
                              color=colors[counter % (len(colors) - 1)], label=rxn_idx_n[str(s_idx)])
            counter += 1
    leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})
    leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()

    a_x_left.set_xlabel("Time/sec")

    a_x_left.set_ylabel("R")
    a_x_right.set_ylabel("T/K")

    rxn_idx_str = "_".join(str(x) for x in reaction_idx)
    plt.title("reaction rates and Temp")

    fig.savefig(os.path.join(file_dir, "output",
                             "reaction_rate_" + rxn_idx_str + "_" + str(end_t) + ".jpg"), dpi=500)
    plt.close()


def plot_reaction_pair_rate_ratio(file_dir, rxn_idx_pair=None, spe_idx_pair=None, tau=10.0, end_t=1.0, tag="M"):
    """
    plot reaction rates ratios given reaction index pair
    """

    colors, markers, _ = get_colors_markers_linestyles()

    _, rxn_idx_n = psri.parse_reaction_and_its_index(file_dir)

    if rxn_idx_pair is None:
        rxn_idx_pair = OrderedDict()
    if spe_idx_pair is None:
        spe_idx_pair = OrderedDict()

    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    rxn_rates = np.loadtxt(os.path.join(file_dir, "output",
                                        "reaction_rate_dlsode_" + str(tag) + ".csv"), delimiter=",")
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_" + str(tag) + ".csv"), delimiter=",")

    # the time point where reference time tau is
    tau_time_point = float(tau) / time[-1] * len(time)
    end_point = int(end_t * tau_time_point)
    if end_point >= len(time):
        end_point = len(time) - 1
    delta_n = int(end_point / 25)
    if delta_n is 0:
        delta_n = 1

    # specify label for lines
    labels = [str(rxn_idx_n[i]) for i in rxn_idx_pair]

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    for idx, (r_idx, s_idx) in enumerate(zip(rxn_idx_pair, spe_idx_pair)):
        data_y = (rxn_rates[::delta_n, int(r_idx)] * conc[::delta_n, int(spe_idx_pair[s_idx])])  \
            / (rxn_rates[::delta_n, rxn_idx_pair[r_idx]] * conc[::delta_n, int(s_idx)])
        a_x.semilogy(time[::delta_n], data_y,
                     color=colors[idx % len(colors)], marker=markers[idx % len(markers)], label=labels[idx % len(labels)])

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 10.0})
    leg.get_frame().set_alpha(0.7)

    a_x.set_xlim([time[0], time[end_point]])
    a_x.grid()

    a_x.set_xlabel("Time/s")
    a_x.set_ylabel("Ratio")
    a_x.set_title("")

    fig.tight_layout()
    # figure name
    fig_name = "reaction_rate_pair_ratios.jpg"
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_spe_path_prob(file_dir, top_n=10, exclude_names=None, init_spe=62, atom_followed="C",
                       end_t=1.0, end_spe=62, species_path=False):
    """
    plot spe_path_prob give species name 
    """
    if exclude_names is None:
        exclude_names = []

    colors, markers, _ = get_colors_markers_linestyles()
    if species_path is False:
        prefix = ""
    else:
        prefix = "species_"
    d_f = pattern_statistics.path_prob_terminating_with_spe(file_dir, init_spe=init_spe,
                                                            atom_followed=atom_followed, end_t=end_t,
                                                            end_spe=end_spe, species_path=species_path)
    data = [float(x) for x in d_f['frequency']][0:top_n]
    data_c = deepcopy(data)
    for idx, _ in enumerate(data_c):
        if idx >= 1:
            data_c[idx] += data_c[idx - 1]
    delta_n = int(len(data_c) / 25)
    if delta_n is 0:
        delta_n = 1
    data_c = data_c[::delta_n]

    s_idx_n, _ = psri.parse_spe_info(file_dir)
    spe_name = s_idx_n[str(end_spe)]

    spe_conc = trajectory.get_normalized_concentration_at_time(
        file_dir, tag="M", end_t=end_t, exclude_names=exclude_names, renormalization=True)
    spe_conc = trajectory.convert_concentration_to_path_prob(
        file_dir, atom_followed=atom_followed, spe_conc=spe_conc, renormalization=True)
    # data_c = trajectory.convert_path_prob_to_concentration(
    #     FILE_DIR, atom_followed=atom_followed, path_prob=data_c)
    spe_conc_const = spe_conc[int(end_spe)]

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    a_x_right = a_x_left.twinx()

    a_x_left.plot(data_c, color=colors[0],
                  marker=markers[0], label="pathway prob")

    a_x_left.plot([0, len(data_c) - 1], [spe_conc_const, spe_conc_const],
                  color=colors[-1], marker=markers[1], label="exact")

    a_x_right.plot((spe_conc_const - data_c) / spe_conc_const,
                   color=colors[1], marker=markers[2], label="Error Percentage")

    leg_left = a_x_left.legend(loc=10, fancybox=True, prop={'size': 10.0})
    leg_right = a_x_right.legend(loc=8, fancybox=True, prop={'size': 10.0})

    a_x_left.yaxis.get_major_formatter().set_powerlimits((0, 1))
    a_x_right.yaxis.get_major_formatter().set_powerlimits((0, 1))
    a_x_left.ticklabel_format(useOffset=False)
    a_x_right.ticklabel_format(useOffset=False)

    x_ticks = [x for x in range(len(data_c))]
    x_tick_labels = [str(delta_n * x + 1) for x in x_ticks]
    a_x_left.set_xticks(x_ticks)
    # a_x_left.set_xticklabels(x_tick_labels, rotation='vertical')
    a_x_left.set_xticklabels(x_tick_labels, rotation=45)

    a_x_left.set_xlim([0 - 0.5, len(data_c) - 0.5])

    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()

    a_x_left.set_xlabel("#path")
    a_x_left.set_ylabel("$\Sigma P_i$ and $\widetilde{X}$")
    a_x_right.set_ylabel("Error")
    ytick_vals = a_x_right.get_yticks()
    a_x_right.set_yticklabels(['{:3.2f}%'.format(x * 100) for x in ytick_vals])

    a_x_left.set_title(spe_name + " @" + str(end_t) + " tau")

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output",
                             prefix + "path_prob_cumulative_" + spe_name + "_" + str(end_t) + ".jpg"), dpi=500)
    plt.close()


def plot_top_n_spe_concentration(file_dir, exclude_names=None, atom_followed="C", end_t=0.5, top_n=10):
    """
    plot_top_n_spe_concentration
    """
    if exclude_names is None:
        exclude_names = []
    spe_conc = trajectory.get_normalized_concentration_at_time(
        file_dir, tag="M", end_t=end_t, exclude_names=exclude_names, renormalization=True)
    spe_conc = trajectory.convert_concentration_to_path_prob(
        file_dir, atom_followed=atom_followed, spe_conc=spe_conc, renormalization=True)

    s_idx_n, _ = psri.parse_spe_info(file_dir)

    spe_top_n_idx_list = spe_conc.argsort()[-top_n:][::-1]

    top_n_spe_conc = [spe_conc[x] for x in spe_top_n_idx_list]
    top_n_spe_name = [s_idx_n[str(x)] for x in spe_top_n_idx_list]

    # figure name
    fig_name = "top_n_spe_concentration_" + \
        str(end_t) + "_" + str(top_n) + ".jpg"

    # specify label for lines
    bins = np.array([x for x in range(len(top_n_spe_conc))])
    height = top_n_spe_conc
    labels = top_n_spe_name

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    rects = a_x.bar(bins, height=height)

    # add number above
    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            a_x.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                     '%.2e' % float(height),
                     ha='center', va='bottom', size=6.0)
    autolabel(rects)

    a_x.set_xticks(bins)
    a_x.set_xticklabels(labels, rotation=-15)

    a_x.get_yaxis().set_visible(False)

    a_x.spines['right'].set_visible(False)
    a_x.spines['top'].set_visible(False)
    a_x.spines['left'].set_visible(False)
    a_x.spines['right'].set_visible(False)

    # a_x.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='off')  # labels along the bottom edge are off

    a_x.set_title("Time = " + str(end_t) + "$\\tau$")

    # background
    a_x.axis('on')
    # a_x.axis('off')

    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.30, 0.30, 0.6, 0.5]
    a_x_2 = fig.add_axes([left, bottom, width, height])

    # specify label for lines
    bins = bins[0:-1]
    height = top_n_spe_conc[1::]
    labels = top_n_spe_name[1::]

    rects_2 = a_x_2.bar(bins, height=height, color='r')

    # add number above
    def autolabel_2(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            a_x_2.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                       '%.2e' % float(height),
                       ha='center', va='bottom', size=6.0, color='r')
    autolabel_2(rects_2)

    a_x_2.axis('off')

    # fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)


def plot_rxn_rate_constant(file_dir):
    """
    plot reaction rate constant, read data from files
    """

    colors, markers, _ = get_colors_markers_linestyles()

    beta = np.loadtxt(os.path.join(file_dir, "output",
                                   "beta.csv"), dtype=float, delimiter=",")
    rxn_name = np.loadtxt(os.path.join(
        file_dir, "output", "rxn_name.csv"), dtype=str, delimiter=",")
    rxn_rate_constant = np.loadtxt(os.path.join(
        file_dir, "output", "rate_constant.csv"), dtype=float, delimiter=",")

    # combine duplicated reactions, same reaction name
    rxn_name_map = dict()
    for idx, val in enumerate(rxn_name):
        if val not in rxn_name_map:
            rxn_name_map[val] = [idx]
        else:
            rxn_name_map[val].append(idx)

    # print(rxn_name_map)
    rxn_name_list = ['O2 + npropyl <=> npropyloo', 'O2 + npropyl <=> C3H6 + HO2',
                     'O2 + npropyl <=> QOOH_1', 'O2 + npropyl <=> OH + propoxide',  'O2 + npropyl <=> QOOH_2']
    # rxn_name_list = ['O2 + ipropyl <=> ipropyloo', 'O2 + ipropyl <=> C3H6 + HO2',
    #                  'O2 + ipropyl <=> OH + propoxide', 'O2 + ipropyl <=> QOOH_3']

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=True)
    counter = 0
    for _, r_n in enumerate(rxn_name_list):
        for idx2, val2 in enumerate(rxn_name_map[r_n]):
            if idx2 == 0:
                y_value = np.array(rxn_rate_constant[:, int(val2)])
            else:
                y_value += np.array(rxn_rate_constant[:, int(val2)])
        counter += 1
        a_x.semilogy(
            beta, y_value, color=colors[counter], marker=markers[counter], label=r_n)

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 10.0})
    leg.get_frame().set_alpha(0.7)
    a_x.grid()

    a_x.set_xlabel("1000/T(K$^{-1}$)")
    a_x.set_ylabel("k(cm$^{3}$ molecule$^{-1}$s$^{-1}$)")
    a_x.set_title("O$_2$ + npropyl")

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output",
                             "reaction_rate_constant.jpg"), dpi=500)
    plt.close()


def plot_pathway_prob_vs_time(file_dir, init_spe=62, atom_followed="C", end_t=1.0,
                              top_n=1, species_path=True):
    """
    plot pathway probability vs. time
    """
    if top_n is None:
        top_n = 1
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    # f_n_pn = os.path.join(file_dir, "output",
    #                       prefix + "pathway_name_selected" + suffix + ".csv")
    f_n_pt = os.path.join(file_dir, "output",
                          prefix + "pathway_time_candidate" + suffix + ".csv")
    f_n_pp = os.path.join(file_dir, "output",
                          prefix + "pathway_prob" + suffix + ".csv")
    # data_pn = np.loadtxt(f_n_pn, dtype=str, delimiter=",")
    data_pt = np.loadtxt(f_n_pt, dtype=float, delimiter=",")
    data_pp = np.loadtxt(f_n_pp, dtype=float, delimiter=",")

    colors, markers, _ = get_colors_markers_linestyles()

    dim_n = len(np.shape(data_pp))

    if dim_n == 1:
        data_x = data_pt
        data_y = data_pp
    elif dim_n == 2:
        data_x = data_pt[0, :]
        data_y = data_pp[0:top_n, :]

    data_x = data_x * end_t

    # labels = data_pn[0:top_n]
    labels = ["P" + str(i + 1) for i in range(top_n)]
    fig_name = prefix + "pathway_prob_vs_time" + suffix + ".jpg"

    delta_n = int(len(data_x) / 50)
    if delta_n is 0:
        delta_n = 1

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    for idx in range(len(data_y)):
        a_x.plot(data_x, data_y[idx, :],
                 color=colors[idx % len(colors)],
                 marker=markers[idx % len(markers)],
                 label=labels[idx % len(labels)],
                 markevery=delta_n)

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 15.0})
    leg.get_frame().set_alpha(0.7)

    y_vals = a_x.get_yticks()
    a_x.set_yticklabels(['{:.1e}'.format(x) for x in y_vals])

    a_x.set_xlim([data_x[0], data_x[-1]])
    a_x.grid()

    a_x.set_xlabel("Time/$\\tau$", fontsize=15.0)
    a_x.set_ylabel("Pathway Probability", fontsize=15.0)
    # a_x.set_title("CO", fontsize=15.0)

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_pathway_AT(file_dir, init_spe=62, atom_followed="C", end_t=1.0,
                    path_idx=0, species_path=True):
    """
    plot pathway arrival time
    """
    if path_idx is None:
        path_idx = 0
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    f_n_pn = os.path.join(file_dir, "output",
                          prefix + "pathway_name_candidate" + suffix + ".csv")
    f_n_pa = os.path.join(file_dir, "output",
                          prefix + "pathway_AT" + suffix + ".csv")
    data_pn = np.loadtxt(f_n_pn, dtype=str, delimiter=",")
    data_pa = np.loadtxt(f_n_pa, dtype=float, delimiter=",")

    fig_name = prefix + "pathway_AT" + \
        suffix + "_" + str(path_idx + 1) + ".jpg"

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    # arguments are passed to np.histogram
    data_hist = data_pa[path_idx, :]
    weights = np.ones_like(data_hist) / float(len(data_hist))
    a_x.hist(data_hist, bins=75, weights=weights, facecolor='green')

    a_x.grid()

    a_x.set_xlabel("Time/s")
    a_x.set_ylabel("Probability")
    a_x.set_title(data_pn[path_idx])

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_pathway_AT_no_IT(file_dir, init_spe=62, atom_followed="C", end_t=1.0,
                          path_idx=0, species_path=True):
    """
    plot pathway arrival time
    """
    if path_idx is None:
        path_idx = 0
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    f_n_pn = os.path.join(file_dir, "output",
                          prefix + "pathway_name_candidate" + suffix + ".csv")
    f_n_pa = os.path.join(file_dir, "output",
                          prefix + "pathway_AT_no_IT" + suffix + ".csv")
    data_pn = np.loadtxt(f_n_pn, dtype=str, delimiter=",")
    data_pa = np.loadtxt(f_n_pa, dtype=float, delimiter=",")

    fig_name = prefix + "pathway_AT_no_IT" + \
        suffix + "_" + str(path_idx + 1) + ".jpg"

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    # arguments are passed to np.histogram
    data_hist = data_pa[path_idx, :]
    weights = np.ones_like(data_hist) / float(len(data_hist))
    a_x.hist(data_hist, bins=75, weights=weights, facecolor='green')

    a_x.grid()

    a_x.set_xlabel("Time/s")
    a_x.set_ylabel("Probability")

    title_n = data_pn[path_idx]

    # load spe and reaction info
    spe_idx_n, _ = psri.parse_spe_info(file_dir)
    _, idx_rxn_n = psri.parse_reaction_and_its_index(file_dir)
    # convert species reaction index to real species and reactions
    title_n = psri.pathname_to_real_spe_reaction(spe_idx_n, idx_rxn_n, title_n)

    a_x.set_title(title_n, fontsize=6)

    # fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_pathway_AT_with_SP(file_dir, init_spe=62, atom_followed="C", end_t=1.0,
                            path_idx=0, species_path=True):
    """
    plot pathway arrival time with terminal species survial probability
    """
    if path_idx is None:
        path_idx = 0
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    f_n_pn = os.path.join(file_dir, "output",
                          prefix + "pathway_name_candidate" + suffix + ".csv")
    f_n_pa = os.path.join(file_dir, "output",
                          prefix + "pathway_AT_with_SP" + suffix + ".csv")
    f_n_psp = os.path.join(file_dir, "output",
                           prefix + "pathway_SP" + suffix + ".csv")
    data_pn = np.loadtxt(f_n_pn, dtype=str, delimiter=",")
    data_pa = np.loadtxt(f_n_pa, dtype=float, delimiter=",")
    data_psp = np.loadtxt(f_n_psp, dtype=float, delimiter=",")

    fig_name = prefix + "pathway_AT_with_SP" + \
        suffix + "_" + str(path_idx + 1) + ".jpg"

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    # arguments are passed to np.histogram
    data_hist = data_pa[path_idx, :]
    data_weights = data_psp[path_idx, :]
    weights = data_weights / float(np.sum(data_weights))
    a_x.hist(data_hist, bins=75, weights=weights, facecolor='green')

    a_x.grid()

    a_x.set_xlabel("Time/s")
    a_x.set_ylabel("Probability")
    a_x.set_title(data_pn[path_idx])

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_first_passage_time(file_dir, init_spe=62, atom_followed="C", end_t=1.0,
                            path_idx=0, species_path=True):
    """
    plot pathway arrival time
    """
    if path_idx is None:
        path_idx = 0
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    f_n_pn = os.path.join(file_dir, "output",
                          prefix + "pathway_name_candidate" + suffix + ".csv")
    f_n_pa = os.path.join(file_dir, "output",
                          prefix + "pathway_AT" + suffix + ".csv")
    data_pn = np.loadtxt(f_n_pn, dtype=str, delimiter=",")
    data_pa = np.loadtxt(f_n_pa, dtype=float, delimiter=",")

    fig_name = prefix + "first_passage_time" + \
        suffix + "_" + str(path_idx + 1) + ".jpg"

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)
    # arguments are passed to np.histogram
    data_hist = data_pa[path_idx, :]
    average_time = np.average(data_hist)

    weights = np.ones_like(data_hist) / float(len(data_hist))
    a_x.hist(data_hist, bins=100, weights=weights, facecolor='blue')

    a_x.grid()

    a_x.set_xlabel("Time/s")
    a_x.set_ylabel("Probability")
    a_x.set_title("S" + re.search(r'\d+', data_pn[path_idx]).group())

    xmin, xmax = a_x.get_xlim()
    ymin, ymax = a_x.get_ylim()
    a_x.text(xmin + 0.1 * (xmax - xmin), ymin + 0.9 * (ymax - ymin),
             "Average Passage Time:  " + '{:1.3f}'.format(float(average_time)))

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_Merchant_f_value(file_dir, init_spe=62, atom_followed="C",
                          begin_t=0.0, end_t=1.0, tau=10.0,
                          species_path=True, spe_idx=None, path_idx=None):
    """
    plot pathway arrival time
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    id_tmp = ""
    if spe_idx is None or spe_idx is []:
        id_tmp = ""
    elif isinstance(spe_idx, int):
        id_tmp = str(spe_idx)
    else:
        for x_t in spe_idx:
            if id_tmp == "":
                id_tmp = str(x_t)
            else:
                id_tmp += "_" + str(x_t)

    s_idx_n, _ = psri.parse_spe_info(file_dir)

    f_n_p_p = os.path.join(file_dir, "output", prefix +
                           "pathway_prob" + suffix + ".csv")
    path_prob = np.loadtxt(f_n_p_p, dtype=float, delimiter=',')

    dim = len(np.shape(path_prob))
    if dim != 2:
        return
    # at least 2 points
    n_points = np.shape(path_prob)[1]
    if n_points < 2:
        return

    if id_tmp != "":
        suffix += "_" + id_tmp

    f_n_s_p_c = os.path.join(
        file_dir, "output", prefix + "pathway_species_production_count" + suffix + ".csv")
    spe_numbers = np.loadtxt(f_n_s_p_c, dtype=float, delimiter=',')

    f_full = [np.dot(path_prob[:, col], spe_numbers)
              for col in range(np.shape(path_prob)[1])]
    print(f_full)

    time_v = np.linspace(begin_t * tau, end_t * tau, n_points + 1)
    time_v = time_v[1::]
    print(time_v)

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    delta_1 = int(len(f_full) / 5)
    a_x.plot(time_v, f_full, label=s_idx_n[str(
        init_spe)], marker='+', markevery=delta_1)

    if path_idx is not None and isinstance(path_idx, list):
        spe_numbers_selected = np.zeros(len(spe_numbers))
        for p_idx in path_idx:
            spe_numbers_selected[p_idx] = spe_numbers[p_idx]
        f_selected = [np.dot(path_prob[:, col], spe_numbers_selected)
                      for col in range(np.shape(path_prob)[1])]
        delta_2 = int(len(f_selected) / 5)
        a_x.plot(time_v, f_selected, label="selected pathways",
                 marker='o', markevery=delta_2)

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 15.0})
    leg.get_frame().set_alpha(0.7)

    y_vals = a_x.get_yticks()
    a_x.set_yticklabels(['{:.2f}'.format(x) for x in y_vals])

    a_x.set_xlim([time_v[0], time_v[-1]])
    a_x.grid()

    a_x.set_xlabel("Time/s", fontsize=15.0)
    a_x.set_ylabel("f", fontsize=15.0)
    a_x.set_title("Merchant f", fontsize=15.0)

    fig.tight_layout()
    if path_idx is not None:
        suffix += "_selected"
    fig_name = prefix + "Merchant_f_vs_time" + suffix + ".jpg"
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_Merchant_alpha_value_vs_time(file_dir, init_spe=10, atom_followed="C", end_t=1.0, species_path=False,
                                      s_idx=10, r_idx=736):
    """
    calculate Merchat alpha value at time point as in time.csv, not at time zero
    """
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    prefix = ""
    if species_path is True:
        prefix = "species_"
    time_v = np.loadtxt(os.path.join(file_dir, "output",
                                     "time_dlsode_M.csv"), dtype=float, delimiter=',')

    merchant_alpha_fn = os.path.join(
        file_dir, "output", prefix + "Merchant_alpha_" + "S" + str(s_idx) + "_R" + str(r_idx) + suffix + ".csv")
    merchant_alpha_v = np.loadtxt(
        merchant_alpha_fn, dtype=float, delimiter=',')

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    a_x.plot(time_v, merchant_alpha_v, label='$\\alpha$')

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 15.0})
    leg.get_frame().set_alpha(0.7)

    y_vals = a_x.get_yticks()
    a_x.set_yticklabels(['{:.2f}'.format(x) for x in y_vals])

    a_x.set_xlim([time_v[0], time_v[-1]])
    a_x.grid()

    a_x.set_xlabel("Time/s", fontsize=15.0)
    a_x.set_ylabel("$\\alpha$", fontsize=15.0)
    a_x.set_title("Merchant $\\alpha$", fontsize=15.0)

    fig.tight_layout()
    fig_name = prefix + "Merchant_alpha_vs_time" + "_S" + \
        str(s_idx) + "_R" + str(r_idx) + suffix + ".jpg"
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


def plot_Merchant_alpha_and_f_value(file_dir, init_spe=62, atom_followed="C",
                                    begin_t=0.0, end_t=1.0, tau=10.0,
                                    species_path=True, s_idx=10, r_idx=736):
    """
    plot pathway arrival time
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"
    suffix = naming.get_suffix(file_dir, init_spe=init_spe,
                               atom_followed=atom_followed, end_t=end_t)
    suffix_ref = deepcopy(suffix)
    id_tmp = str(s_idx)

    s_idx_n, _ = psri.parse_spe_info(file_dir)

    f_n_p_p = os.path.join(file_dir, "output", prefix +
                           "pathway_prob" + suffix + ".csv")
    p_1 = np.loadtxt(f_n_p_p, dtype=float, delimiter=',')

    dim = len(np.shape(p_1))
    if dim != 2:
        return
    # at least 2 points
    n_points = np.shape(p_1)[1]
    if n_points < 2:
        return

    if id_tmp != "":
        suffix += "_" + id_tmp

    f_n_s_p_c = os.path.join(
        file_dir, "output", prefix + "pathway_species_production_count" + suffix + ".csv")
    p_2 = np.loadtxt(f_n_s_p_c, dtype=float, delimiter=',')

    merchant_f_v = [np.dot(p_1[:, col], p_2)
                    for col in range(np.shape(p_1)[1])]
    print(merchant_f_v)

    time_v = np.linspace(begin_t * tau, end_t * tau, n_points + 1)
    time_v = time_v[1::]
    print(time_v)

    time_ref = np.loadtxt(os.path.join(file_dir, "output",
                                       "time_dlsode_M.csv"), dtype=float, delimiter=',')
    merchant_alpha_fn = os.path.join(
        file_dir, "output", prefix + "Merchant_alpha_" + "S" + str(s_idx) + "_R" + str(r_idx) + suffix_ref + ".csv")
    alpha_ref = np.loadtxt(
        merchant_alpha_fn, dtype=float, delimiter=',')

    alpha_f_v = np.zeros(n_points)
    for i in range(n_points):
        alpha_f_v[i] = merchant_f_v[i] * \
            interpolation.interp1d(time_ref, alpha_ref, time_v[i])

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=False)

    delta_n1 = int(len(time_v) / 25)
    a_x.plot(time_v, merchant_f_v, label='$f$', marker='*', markevery=delta_n1)
    delta_n2 = int(len(time_ref) / 25)
    a_x.plot(time_ref, alpha_ref, label='$\\alpha$',
             markevery=delta_n2, marker='o')
    a_x.plot(time_v, alpha_f_v, label='$\\alpha$ corrected $f$',
             marker='+', markevery=delta_n1)

    leg = a_x.legend(loc=0, fancybox=True, prop={'size': 15.0})
    leg.get_frame().set_alpha(0.7)

    y_vals = a_x.get_yticks()
    a_x.set_yticklabels(['{:.2f}'.format(x) for x in y_vals])

    a_x.set_xlim([time_v[0], time_v[-1]])
    a_x.grid()

    a_x.set_xlabel("Time/s", fontsize=15.0)
    a_x.set_ylabel("f", fontsize=15.0)
    a_x.set_title("Merchant f", fontsize=15.0)

    fig.tight_layout()
    fig_name = prefix + "Merchant_alpha_and_f_vs_time" + suffix + ".jpg"
    fig.savefig(os.path.join(file_dir, "output", fig_name), dpi=500)
    plt.close()


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    G_S = global_settings.get_setting(FILE_DIR)

    SPE_LIST = [59, 10, 12, 13, 59]
    # SPE_LIST, _, _ = trajectory.get_species_with_top_n_concentration(
    #     FILE_DIR, exclude=None, top_n=10,
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], end_t=G_S['end_t'],
    #     tag="M", atoms=["C", "O"])
    plot_concentrations(FILE_DIR, spe_idx=SPE_LIST,
                        tau=G_S['tau'], end_t=0.01, tag="M", exclude_names=None,
                        renormalization=False, semilogy=True, hasTemp=True)
    # plot_concentrations(FILE_DIR, spe_idx=SPE_LIST,
    #                     tau=G_S['tau'], end_t=1.0, tag="M", exclude_names=None,
    #                     renormalization=False, semilogy=True, hasTemp=True)
    # plot_concentrations(FILE_DIR, spe_idx=[10],
    #                     tau=G_S['tau'], end_t=1.0, tag="M", exclude_names=None,
    #                     renormalization=False, semilogy=True, hasTemp=True)
    # plot_spe_concentrations_derivative(FILE_DIR, spe_idx=[62, 14, 15, 59],
    #                                    tau=G_S['tau'], end_t=0.95, tag="M",
    #                                    exclude_names=None, renormalization=False)
    # SPE_LIST = [14, 59, 17, 44, 38, 86,  69, 15, 82]
    # for es in SPE_LIST:
    #     plot_path_length_statistics(
    #         FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=0.9, end_spe=es)
    # plot_pathway_prob_vs_time(
    #     FILE_DIR, init_spe=60, atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #     top_n=6, species_path=G_S['species_path'])
    # plot_spe_drc(FILE_DIR, spe_idx=[25, 39, 45, 60, 61, 72, 54],
    #              tau=G_S['tau'], end_t=0.80, tag="M", reciprocal=True)
    # plot_chattering_group_drc(
    #     FILE_DIR, tau=G_S['tau'], end_t=0.80, tag="M", reciprocal=True)
    # plot_spe_drc(FILE_DIR, spe_idx=[62, 94, 101, 46, 17, 16, 14, 44, 47],
    #              tau=G_S['tau'], end_t=1.0, tag="M", reciprocal=True)
    # plot_spe_drc(FILE_DIR, spe_idx=[60, 61, 87, 90, 94, 101],
    #              tau=G_S['tau'], end_t=0.25, tag="M", reciprocal=True)
    # plot_chattering_group_drc(
    #     FILE_DIR, tau=G_S['tau'], end_t=0.25, tag="M", reciprocal=True)
    # for p_i in range(10):
    #     plot_pathway_AT(
    #         FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #         path_idx=p_i, species_path=True)
    # for p_i in range(20):
    #     plot_pathway_AT_no_IT(
    #         FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #         path_idx=p_i, species_path=True)
    # for p_i in range(20):
    #     plot_pathway_AT_with_SP(
    #         FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #         path_idx=p_i, species_path=True)
    # for p_i in range(len(G_S['end_s_idx'])):
    #     plot_first_passage_time(
    #         FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #         path_idx=p_i, species_path=True)

    # PATH_IDX = [0, 1, 2, 3, 4, 6, 11, 44, 59, 66,
    #             68, 93, 115, 138, 153, 165, 166, 245, 477]
    # plot_Merchant_f_value(FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                       begin_t=G_S['begin_t'], end_t=G_S['end_t'], tau=G_S['tau'],
    #                       species_path=G_S['species_path'], spe_idx=[10],
    #                       path_idx=PATH_IDX)

    # TIME_AXIS = pathway_time_2_array_index(
    #     FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #     species_path=G_S['species_path'], time=G_S['mc_t'])
    # plot_cumulative_pathway_prob(FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                              end_t=G_S['end_t'], tau=G_S['tau'], species_path=G_S['species_path'],
    #                              top_n=25, time_axis=TIME_AXIS)

    # plot_Merchant_alpha_value_vs_time(
    #     FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'], end_t=G_S['end_t'],
    #     species_path=G_S['species_path'], s_idx=10, r_idx=736)

    # plot_Merchant_alpha_and_f_value(FILE_DIR, init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                                 begin_t=G_S['begin_t'], end_t=G_S['end_t'], tau=G_S['tau'],
    #                                 species_path=G_S['species_path'], s_idx=10)
