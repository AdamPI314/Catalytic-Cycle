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


def plot_concentrations(file_dir, spe_idx=None):
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

    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    # linestyles = Line2D.lineStyles.keys()

    s_n_idx = {"-1": "Temp", "0": "HE", "10": "OH", "12": "HO2",
               "62": "C3H8", "9": "O2", "2": "N2", "1": "Ar",
               "13": "H2O2", "17": "CH2O",
               "46": "CH2CHO", "50": "CH3CH2OO",
               "51": "CH3CH2OOH", "52": "CH2CH2OOH"}
    if spe_idx is None:
        spe_idx = [0]
    time = np.loadtxt(os.path.join(
        file_dir, "output", "time_dlsode_M.csv"), delimiter=",")
    conc = np.loadtxt(os.path.join(file_dir, "output",
                                   "concentration_dlsode_M.csv"), delimiter=",")
    temp = np.loadtxt(os.path.join(file_dir, "output",
                                   "temperature_dlsode_M.csv"), delimiter=",")

    counter = 0

    fig, a_x_left = plt.subplots(1, 1, sharex=True, sharey=False)
    for s_idx in spe_idx:
        if s_idx == -1:
            a_x_right = a_x_left.twinx()
            a_x_right.plot(time, temp,
                           color=colors[counter % len(colors)], label=s_n_idx[str(s_idx)])
        else:
            a_x_left.semilogy(time, conc[:, s_idx],
                              color=colors[counter % len(colors)], label=s_n_idx[str(s_idx)])
        counter += 1
    a_x_left.legend(loc=0, prop={'size': 10.0})
    a_x_right.legend(loc=0, prop={'size': 10.0})

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
    # plot_concentrations(FILE_DIR, [13, 17, 62, -1])
    # plot_concentrations(FILE_DIR, [46, 50, 51, 52, -1])
    print(FILE_DIR)
