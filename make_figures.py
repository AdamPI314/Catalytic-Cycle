"""
plot jobs
"""
import os
import sys
from copy import deepcopy
import numpy as np
from collections import OrderedDict

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
    end_point = int(tau * len(time))
    delta_n = int(end_point / 25)
    if delta_n is 0:
        delta_n = 1

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
    end_point = int(tau * len(time))
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
                             "reaction_rate_" + rxn_idx_str + "_" + str(tau) + ".jpg"), dpi=500)
    plt.close()


def plot_reaction_pair_rate_ratio(file_dir, rxn_idx_pair=None, spe_idx_pair=None, tau=1.0, tag="M"):
    """
    plot reaction rates ratios given reaction index pair
    """

    colors, markers, _ = get_colors_markers_linestyles()

    _, rxn_idx_n = psri.parse_reaction_and_its_index(os.path.join(
        file_dir, "output", "reaction_labelling.csv"))

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

    end_point = int(tau * len(time))
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
                     color=colors[idx], marker=markers[idx], label=labels[idx])

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


def plot_spe_path_prob(file_dir, spe_name="C3H8", top_n=10, exclude_names=None, init_spe=62, atom_followed="C", tau=1.0, pathwayEndWith="ALL"):
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
    print(spe_name, "#path:\t", len(data_c))
    delta_n = int(len(data_c) / 25)
    if delta_n is 0:
        delta_n = 1
    data_c = data_c[::delta_n]

    _, s_n_idx = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))

    spe_conc = trajectory.get_normalized_concentration_at_time(
        file_dir, tag="M", tau=tau, exclude_names=exclude_names, renormalization=True)
    spe_conc = trajectory.convert_concentration_to_path_prob(
        file_dir, atom_followed=atom_followed, spe_conc=spe_conc, renormalization=True)
    # data_c = trajectory.convert_path_prob_to_concentration(
    #     FILE_DIR, atom_followed=atom_followed, path_prob=data_c)
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


def plot_top_n_spe_concentration(file_dir, exclude_names=None, atom_followed="C", tau=0.5, top_n=10):
    """
    plot_top_n_spe_concentration
    """
    if exclude_names is None:
        exclude_names = []
    spe_conc = trajectory.get_normalized_concentration_at_time(
        file_dir, tag="M", tau=tau, exclude_names=exclude_names, renormalization=True)
    spe_conc = trajectory.convert_concentration_to_path_prob(
        file_dir, atom_followed=atom_followed, spe_conc=spe_conc, renormalization=True)

    s_idx_n, _ = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))

    spe_top_n_idx_list = spe_conc.argsort()[-top_n:][::-1]

    top_n_spe_conc = [spe_conc[x] for x in spe_top_n_idx_list]
    top_n_spe_name = [s_idx_n[str(x)] for x in spe_top_n_idx_list]

    # figure name
    fig_name = "top_n_spe_concentration_" + \
        str(tau) + "_" + str(top_n) + ".jpg"

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

    a_x.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')  # labels along the bottom edge are off

    a_x.set_title("$\\tau$ = " + str(tau))

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


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # plot_pathway_prob(FILE_DIR, tau=0.2)
    G_S = global_settings.get_setting()

    SPE_IDX, SPE_NAMES, SPE_EXCLUDE_NAME = trajectory.get_species_with_top_n_concentration(
        FILE_DIR, exclude=None, top_n=G_S['top_n_s'], tau=G_S['tau'], tag=G_S['tag'], atoms=[G_S['atom_f']])
    # plot_concentrations(
    #     FILE_DIR, spe_idx=[62, 17, 66, 86], tag=G_S['tag'],
    #     exclude_names=SPE_EXCLUDE_NAME, renormalization=True)
    # plot_concentrations(
    #     FILE_DIR, spe_idx=SPE_IDX, tag=G_S['tag'],
    #     exclude_names=SPE_EXCLUDE_NAME, renormalization=True)
    # plot_reaction_rates(
    #     FILE_DIR, reaction_idx=[1068, 1070, 1072, 1074, 1076], tau=1.0, tag=G_S['tag'])
    # for spe_n in SPE_NAMES:
    #     plot_spe_path_prob(FILE_DIR, spe_name=spe_n, top_n=G_S['top_n_p'],
    #                        exclude_names=SPE_EXCLUDE_NAME, tau=G_S['tau'])
    # plot_rxn_rate_constant(FILE_DIR)
    # R_IDX_PAIR, S_IDX_PAIR = global_settings.get_fast_rxn_trapped_spe(FILE_DIR)
    # plot_reaction_pair_rate_ratio(
    #     FILE_DIR, rxn_idx_pair=R_IDX_PAIR, spe_idx_pair=S_IDX_PAIR, tau=1.0, tag="M")
    # plot_top_n_spe_concentration(
    #     FILE_DIR, exclude_names=None, atom_followed=G_S['atom_f'], tau=G_S['tau'], top_n=10)
