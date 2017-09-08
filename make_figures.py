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
from matplotlib.lines import Line2D

import parse_spe_reaction_info as psri
import trajectory
import parse_pattern


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

    markers_tmp = []
    for m_k in Line2D.markers:
        try:
            if len(m_k) == 1 and m_k != ' ':
                markers_tmp.append(m_k)
        except TypeError:
            pass
    markers = markers_tmp[2::]
    markers.append(markers_tmp[0])
    markers.append(markers_tmp[1])

    # styles = markers + [
    #     r'$\lambda$',
    #     r'$\bowtie$',
    #     r'$\circlearrowleft$',
    #     r'$\clubsuit$',
    #     r'$\checkmark$']

    colors = ('b', 'g', 'k', 'c', 'm', 'y', 'r')
    # linestyles = Line2D.lineStyles.keys()

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
    markers_tmp = []
    for m_k in Line2D.markers:
        try:
            if len(m_k) == 1 and m_k != ' ':
                markers_tmp.append(m_k)
        except TypeError:
            pass
    markers = markers_tmp[2::]
    markers.append(markers_tmp[0])
    markers.append(markers_tmp[1])

    # styles = markers + [
    #     r'$\lambda$',
    #     r'$\bowtie$',
    #     r'$\circlearrowleft$',
    #     r'$\clubsuit$',
    #     r'$\checkmark$']

    colors = ('b', 'g', 'k', 'c', 'm', 'y', 'r')
    # linestyles = Line2D.lineStyles.keys()

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
                             "reaction_rate_" + rxn_idx_str + ".jpg"), dpi=500)
    plt.close()


def plot_spe_path_prob(file_dir, spe_name="C3H8", top_n=10, exclude_names=None, tau=1.0, renormalization=True):
    """
    plot spe_path_prob give species name 
    """
    if exclude_names is None:
        exclude_names = []

    markers_tmp = []
    for m_k in Line2D.markers:
        try:
            if len(m_k) == 1 and m_k != ' ':
                markers_tmp.append(m_k)
        except TypeError:
            pass

    markers_tmp = markers_tmp + [
        r'$\lambda$',
        r'$\bowtie$',
        r'$\circlearrowleft$',
        r'$\clubsuit$',
        r'$\checkmark$']

    markers = markers_tmp[2::]
    markers.append(markers_tmp[0])
    markers.append(markers_tmp[1])

    colors = ('b', 'g', 'k', 'c', 'm', 'y', 'r')
    # linestyles = Line2D.lineStyles.keys()

    d_f = parse_pattern.parse_path_prob_terminating_with_spe(
        file_dir, spe_name=spe_name)
    data = [float(x) for x in d_f['frequency']][0:top_n]
    data_c = deepcopy(data)
    for idx, _ in enumerate(data_c):
        if idx >= 1:
            data_c[idx] += data_c[idx - 1]

    _, s_n_idx = psri.parse_spe_info(os.path.join(
        file_dir, "output", "species_labelling.csv"))

    spe_conc = trajectory.get_normalized_concentration(
        file_dir, tag="fraction", exclude_names=exclude_names, renormalization=renormalization)
    spe_conc_const = spe_conc[int(tau * len(spe_conc)), int(s_n_idx[spe_name])]

    fig, a_x = plt.subplots(1, 1, sharex=True, sharey=True)
    a_x.plot(data_c, color=colors[0], marker=markers[0], label="pathway prob")

    a_x.plot([0, len(data_c) - 1], [spe_conc_const, spe_conc_const],
             color=colors[-1], marker=markers[1], label="exact")

    a_x.grid()
    a_x.legend(loc=0, fancybox=True, prop={'size': 10.0})

    fig.savefig(os.path.join(file_dir, "output",
                             "path_prob_cumulative_" + spe_name + ".jpg"), dpi=500)
    plt.close()


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    # plot_pathway_prob(FILE_DIR, tau=0.2)
    # plot_concentrations(FILE_DIR, [62, 17, 66, 86, -1])
    SPE_IDX, SPE_EXCLUDE_NAME = trajectory.get_species_with_top_n_concentration(
        FILE_DIR, exclude=None, top_n=10, tau=0.9, tag="M", atoms=["C"])
    # plot_concentrations(
    #     FILE_DIR, spe_idx=SPE_IDX, tag="M", renormalization=False)
    plot_concentrations(
        FILE_DIR, spe_idx=SPE_IDX, tag="fraction", exclude_names=SPE_EXCLUDE_NAME, renormalization=True)
    # plot_reaction_rates(
    #     FILE_DIR, reaction_idx=[1068, 1070, 1072, 1074, 1076], tag="M")
    # plot_spe_path_prob(FILE_DIR, spe_name="C3H8", top_n=1000,
    #                    exclude_names=SPE_EXCLUDE_NAME, tau=0.9, renormalization=True)
    print(FILE_DIR)
