"""
plot jobs
"""
import os
import sys
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.lines import Line2D
from copy import deepcopy


def plot_pathway_prob(file_dir, max_tau=1.0):
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
    end_point = int(max_tau * len(data))
    a_x.plot(data[1:end_point:1], "-.")
    a_x.plot(data_c[1:end_point:1], "-*")
    fig.savefig(os.path.join(file_dir, "output",
                             "pathway_probability.jpg"), dpi=500)

    return


def plot_concentrations(file_dir, spe_idx=None, max_tau=1.0):
    """
    plot concentrations give species index list
    """
    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                markers.append(m)
        except TypeError:
            pass
    # styles = markers + [
    #     r'$\lambda$',
    #     r'$\bowtie$',
    #     r'$\circlearrowleft$',
    #     r'$\clubsuit$',
    #     r'$\checkmark$']

    colors = ('b', 'g', 'k', 'c', 'm', 'y', 'r')
    # linestyles = Line2D.lineStyles.keys()

    s_n_idx = {"-1": "Temp", "0": "HE", "10": "OH", "12": "HO2",
               "62": "C3H8", "9": "O2", "2": "N2", "1": "Ar",
               "13": "H2O2", "17": "CH2O",
               "46": "CH2CHO", "50": "CH3CH2OO",
               "51": "CH3CH2OOH", "52": "CH2CH2OOH",
               "14": "CO", "15": "CO2"}
    if spe_idx is None:
        spe_idx = [0]
    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_fraction.csv"), delimiter=",")
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_fraction.csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_fraction.csv"), delimiter=",")

    counter = 0
    end_point = int(max_tau * len(time))

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time[0:end_point], temp[0:end_point],
                           color=colors[-1], label=s_n_idx[str(s_idx)])
        else:
            a_x_left.semilogy(time[0:end_point], conc[0:end_point, s_idx],
                              color=colors[counter % len(colors)], label=s_n_idx[str(s_idx)])
            counter += 1
    a_x_left.legend(loc=0, prop={'size': 10.0})
    a_x_right.legend(loc=1, prop={'size': 10.0})

    plt.grid(True)
    plt.xlabel("Time/sec")

    a_x_left.set_ylabel("Log(X)")
    a_x_right.set_ylabel("T/K")

    s_n_str = "_".join(s_n_idx[str(x)] for x in spe_idx)
    plt.title(s_n_str)

    fig.savefig(os.path.join(file_dir, "output",
                             "trajectory_" + s_n_str + ".jpg"), dpi=500)
    plt.close()


if __name__ == '__main__':
    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    plot_concentrations(FILE_DIR, [10, 12, -1])
    plot_concentrations(FILE_DIR, [10, 12, -1])
    plot_concentrations(FILE_DIR, [14, 15, -1])
    plot_concentrations(FILE_DIR, [13, 17, 62, -1])
    plot_concentrations(FILE_DIR, [46, 50, 51, 52, -1])
    plot_pathway_prob(FILE_DIR, max_tau=0.2)
    print(FILE_DIR)
