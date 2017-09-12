"""
plot jobs
"""
import os
import sys
from copy import deepcopy
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt

import parse_spe_reaction_info as psri
import trajectory
import parse_pattern
import global_settings
from tools import get_colors_markers_linestyles


def plot_pathway_prob(file_dir, tau=1.0):
    """
    plot pathway prob
    """
    data = np.genfromtxt(os.path.join(file_dir, "output",
                                      "pathway_prob.csv"), delimiter=",")
    data_c = deepcopy(data)
    data_c[0] = 0.0
    for idx, _ in enumerate(data_c):
        if idx >= 1:
            data_c[idx] += data_c[idx - 1]
    print(data_c[-1])
    # np.savetxt(os.path.join(file_dir, "output",
    #                         "pathway_c_prob.csv"), data_c, fmt="%.15f", delimiter=',')

    fig, a_x = plt.subplots(1, 1, sharex=False, sharey=False)
    end_point = int(tau * len(data))
    a_x.plot(data[1:end_point:1], "-.")
    a_x.plot(data_c[1:end_point:1], "-*")
    a_x.grid()
    fig.savefig(os.path.join(file_dir, "output",
                             "pathway_probability.jpg"), dpi=500)

    return


def plot_concentrations(file_dir, spe_idx=None, tau=1.0, tag="fraction", exclude_names=None, renormalization=True):
    """
    plot concentrations give species index list, if exclude is not None, means we are going
    to renormalize the molelar fraction
    """
    if exclude_names is None:
        exclude_names = []

    colors, markers, _ = get_colors_markers_linestyles()

    s_idx_n, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))
    s_idx_n["-1"] = "Temp"
    spe_idx.append(-1)

    if spe_idx is None:
        spe_idx = [0]
    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_" + str(tag) + ".csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_" + str(tag) + ".csv"), delimiter=",")

    conc = trajectory.get_normalized_concentration(
        file_dir, tag=tag, exclude_names=exclude_names, renormalization=renormalization)

    counter = 0
    delta_n = 20
    end_point = int(tau * len(time))

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point:delta_n], temp[0:end_point:delta_n],
                           color=colors[-1], label=s_idx_n[str(s_idx)])
        else:
            if counter < len(colors) - 1:
                m_k = None
            else:
                m_k = markers[(counter + 1 - len(colors)) % (len(markers))]
            a_x_left.semilogy(time[0:end_point:delta_n], conc[0:end_point:delta_n, s_idx], marker=m_k,
                              color=colors[counter % (len(colors) - 1)], label=s_idx_n[str(s_idx)])
            counter += 1
    leg_left = a_x_left.legend(loc=8, fancybox=True, prop={'size': 10.0})
    leg_right = a_x_right.legend(loc=2, fancybox=True, prop={'size': 10.0})
    leg_left.get_frame().set_alpha(0.7)
    leg_right.get_frame().set_alpha(0.7)
    a_x_left.grid()

    a_x_left.set_xlabel("Time/sec")

    a_x_left.set_ylabel("[X]")
    a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_idx_n[str(x)] for x in spe_idx)
    # plt.title(s_n_str)

    fig.savefig(os.path.join(file_dir, "output",
                             "trajectory_" + s_n_str + ".jpg"), dpi=500)
    plt.close()


def plot_reaction_rates(file_dir, reaction_idx=None, tau=1.0, tag="fraction"):
    """
    plot reaction rates give reaction index list
    """

    colors, markers, _ = get_colors_markers_linestyles()

    _, rxn_idx_n = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))
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
    delta_n = 20
    end_point = int(tau * len(time))

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
                             "reaction_rate_" + rxn_idx_str + "_" + str(tau) + ".jpg"), dpi=500)
    plt.close()


def plot_spe_path_prob(file_dir, spe_name="C3H8", top_n=10, exclude_names=None, init_spe=62, atom_followed="C", tau=1.0, pathwayEndWith="ALL", renormalization=True):
    """
    plot spe_path_prob give species name 
    """
    if exclude_names is None:
        exclude_names = []

    colors, markers, _ = get_colors_markers_linestyles()

    d_f = parse_pattern.parse_path_prob_terminating_with_spe(file_dir, spe_name=spe_name, init_spe=init_spe,
                                                             atom_followed=atom_followed, tau=tau, pathwayEndWith=pathwayEndWith)
    data = [float(x) for x in d_f['frequency']][0:top_n]
    data_c = deepcopy(data)
    for idx, _ in enumerate(data_c):
        if idx >= 1:
            data_c[idx] += data_c[idx - 1]
    print(spe_name, "#path:\t", len(data))
    delta_n = int(len(data_c) / 25)
    if delta_n is 0:
        delta_n = 1
    data_c = data_c[::delta_n]

    _, s_n_idx = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))

    spe_conc = trajectory.get_normalized_concentration_at_time(
        FILE_DIR, tag="M", tau=tau, exclude_names=exclude_names, renormalization=True)
    trajectory.convert_concentration_to_path_prob(
        FILE_DIR, atom_followed=atom_followed, spe_conc=spe_conc, renormalization=True)
    spe_conc_const = spe_conc[int(s_n_idx[spe_name])]

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

    a_x_left.set_title(spe_name + " @" + str(tau) + " tau")

    fig.tight_layout()
    fig.savefig(os.path.join(file_dir, "output",
                             "path_prob_cumulative_" + spe_name + "_" + str(tau) + ".jpg"), dpi=500)
    plt.close()


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


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # plot_pathway_prob(FILE_DIR, tau=0.2)
    # plot_concentrations(FILE_DIR, [62, 17, 66, 86, -1])

    G_S = global_settings.get_setting()

    SPE_IDX, SPE_NAMES, SPE_EXCLUDE_NAME = trajectory.get_species_with_top_n_concentration(
        FILE_DIR, exclude=None, top_n=G_S['top_n_s'], tau=G_S['tau'], tag=G_S['tag'], atoms=[G_S['atom_f']])
    # plot_concentrations(
    # FILE_DIR, spe_idx=SPE_IDX, tag=G_S['tag'],
    # exclude_names=SPE_EXCLUDE_NAME, renormalization=True)
    # plot_reaction_rates(
    # FILE_DIR, reaction_idx=[1068, 1070, 1072, 1074, 1076], tau=G_S['tau'], tag=G_S['tag])
    for spe_n in SPE_NAMES:
        plot_spe_path_prob(FILE_DIR, spe_name=spe_n, top_n=G_S['top_n_p'],
                           exclude_names=SPE_EXCLUDE_NAME, tau=G_S['tau'], renormalization=True)
    # plot_rxn_rate_constant(FILE_DIR)
