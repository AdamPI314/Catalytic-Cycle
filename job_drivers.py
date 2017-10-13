"""
Job drivers
"""

import subprocess
import os
import update_settings as us
import parse_spe_reaction_info as psri
import prepare_path_name_time as ppnt
import parse_pattern as pp
import naming
import trajectory
import make_figures as mf
import global_settings


def update_trapped_species_fast_reaction_setting(file_dir):
    """
    update settings.json, primarily for trapped species and fast reactions
    """
    us.update_trapped_species_fast_reaction_setting(file_dir)


def copy_sohr_files(file_dir):
    """
    copy SOHR files from C++ routine
    """
    naming.copy_sohr_files(file_dir)


def species_count(file_dir, top_n=50, norm=False):
    """
    count species occurence in pathway, multiply by accurate pathwap probability
    """
    pp.species_count(file_dir, top_n=top_n, norm=norm)


def reaction_count(file_dir, top_n=50, norm=False):
    """
    count reaction occurence in pathway, multiply by accurate pathwap probability
    """
    pp.reaction_count(file_dir, top_n=top_n, norm=norm)


def initiation_reaction_count(file_dir, top_n=50, norm=False):
    """
    count initiation reaction occurence in pathway, multiply by accurate pathwap probability
    """
    pp.initiation_reaction_count(file_dir, top_n=top_n, norm=norm)


def species_cycle(file_dir, top_n=50, norm=False):
    """
    count species cycle in pathway, multiply by accurate pathwap probability
    """
    pp.species_cycle(file_dir, top_n=top_n, norm=norm)


def species_production_path(file_dir, spe='OH', top_n=50, norm=False):
    """
    count species production pathway or sub-pathway in pathway,
    multiply by accurate pathwap probability
    """
    pp.species_production_path(file_dir, spe=spe, top_n=top_n, norm=norm)


def species_production_reaction(file_dir, spe='OH', top_n=50, norm=False):
    """
    count species production reactions in pathway,
    multiply by accurate pathwap probability
    """
    pp.species_production_reaction(file_dir, spe=spe, top_n=top_n, norm=norm)


def symbolic_path_2_real_path(file_dir, top_n=50, flag="", end_s_idx=None):
    """
    convert symbolic pathway to real pathway with real species name and real reaction name
    flag indicates a specific job, for example, pathway end time = 1.0, the j-th run,
    any unique symbol shall work
    """
    if flag == "":
        out_file_name = "pathname_prob.csv"
    else:
        out_file_name = "pathname_prob_" + str(flag) + ".csv"
    psri.symbolic_path_2_real_path(
        os.path.join(file_dir, "output", "species_labelling.csv"),
        os.path.join(
            file_dir, "output", "reaction_labelling.csv"),
        os.path.join(
            file_dir, "output", "pathway_stat.csv"),
        os.path.join(
            file_dir, "output", out_file_name),
        top_n, end_s_idx)


def delete_non_dlsode_files(file_dir):
    """
    delete none dlsode files
    """
    os.chdir(file_dir)
    cmd = ["find", "./output", "-type", "f",
           "!", "-name", "*dlsode*", "-delete"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


def make_run(file_dir):
    """
    make run
    """
    os.chdir(file_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


def run_dlsode(file_dir, time):
    """
    Run dlsode
    """
    os.chdir(file_dir)
    us.update_dlsode_setting(file_dir, time)
    make_run(file_dir)


def spe_concentration_at_time_w2f(file_dir, max_tau=10.0, tau=1.0):
    """
    write species concentration at a time to file
    """
    os.chdir(file_dir)
    us.update_spe_concentration_at_time_w2f(file_dir, max_tau=max_tau, tau=tau)
    make_run(file_dir)


def run_mc_trajectory(file_dir, n_traj=1000000, atom_followed="C", init_spe=114, max_tau=10.0, tau=1.0):
    """
    Run mc trajectory
    """
    os.chdir(file_dir)
    us.update_mc_trajectory_setting(
        file_dir, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe, max_tau=max_tau, tau=tau)
    make_run(file_dir)


def evaluate_pathway_probability(file_dir, top_n=5, num_t=1, flag="", n_traj=10000,
                                 atom_followed="C", init_spe=114, traj_end_time=100.0, max_tau=10.0, tau=1.0, top_n_s=10, spe_oriented=True, end_s_idx=None):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(file_dir)

    if spe_oriented is True:
        us.update_eval_path_integral(
            file_dir, top_n=top_n * top_n_s, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe, max_tau=max_tau, tau=tau)
        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                file_dir, exclude=None, top_n=top_n_s, traj_end_time=traj_end_time, max_tau=max_tau, tau=tau, tag="M", atoms=[atom_followed])
        ppnt.prepare_pathway_name(
            file_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx)
        ppnt.prepare_pathway_time(
            file_dir, top_n=top_n * top_n_s, num=num_t, flag=flag, tau=tau)
    else:
        us.update_eval_path_integral(
            file_dir, top_n=top_n, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe, max_tau=max_tau, tau=tau)
        ppnt.prepare_pathway_name(
            file_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx)
        ppnt.prepare_pathway_time(
            file_dir, top_n=top_n, num=num_t, flag=flag, tau=tau)

    make_run(file_dir)


# http://stackoverflow.com/questions/3000724/running-matlab-in-the-background
def make_a_figure(file_dir, ind):
    """
    make a figure
    """
    os.chdir(file_dir)
    matlab_script_dir = os.path.join(
        file_dir, "tools/data_analysis/H2_O2_reaction_network_ODE_solver")
    matlab_script_filename = "plot_concentration_v2"
    matlab_cmd = "cd " + matlab_script_dir + "; " + matlab_script_filename + "(" + str(
        ind) + ")" + "; cd ../../..; exit;"
    cmd = ["nohup", "matlab", "-nosplash", "-nodisplay", "-r", matlab_cmd]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()


def make_figures(file_dir):
    """
    make figures
    """
    for i in range(1, 9):
        make_a_figure(file_dir, i)


def propane_make_figures(file_dir):
    """
    make figures for propane system
    """
    g_s = global_settings.get_setting()

    spe_idx, spe_names, spe_exclude_name = trajectory.get_species_with_top_n_concentration(
        file_dir, exclude=None, top_n=g_s['top_n_s'], traj_end_time=g_s['end_t'],
        max_tau=g_s['max_tau'], tau=g_s['tau'], tag=g_s['tag'], atoms=[g_s['atom_f']])
    mf.plot_concentrations(
        file_dir, spe_idx=spe_idx, max_tau=g_s['max_tau'], tau=g_s['tau'], tag=g_s['tag'],
        exclude_names=spe_exclude_name, renormalization=True)
    mf.plot_reaction_rates(
        file_dir, reaction_idx=[1068, 1070, 1072, 1074, 1076], max_tau=g_s['max_tau'], tau=1.0, tag=g_s['tag'])
    for spe_n in spe_names:
        mf.plot_spe_path_prob(file_dir, spe_name=spe_n, top_n=g_s['top_n_p'],
                              exclude_names=spe_exclude_name, tau=g_s['tau'])
    mf.plot_rxn_rate_constant(file_dir)
    r_idx_pair, s_idx_pair = global_settings.get_fast_rxn_trapped_spe(file_dir)
    mf.plot_reaction_pair_rate_ratio(
        file_dir, rxn_idx_pair=r_idx_pair, spe_idx_pair=s_idx_pair, max_tau=g_s['max_tau'], tau=1.0, tag="M")
    mf.plot_top_n_spe_concentration(
        file_dir, exclude_names=None, atom_followed=g_s['atom_f'], tau=g_s['tau'], top_n=10)


def send_email(file_dir):
    """
    send email to elliot.srbai@gmail.com
    """
    os.chdir(file_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
           "-u", "RUNNING JOB", "-m", "JOB FINISHED." + "\n" + file_dir,
           "-a", "./log.txt"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'ab+')

    try:
        result = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=error_file)
    except subprocess.CalledProcessError as error:
        print(error)
        exit(1)

    if result.stdout is not None:
        out = result.stdout.read()
        out_file.write(out)

    out_file.close()
    error_file.close()
