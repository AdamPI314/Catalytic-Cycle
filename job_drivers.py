"""
Job drivers
"""

import subprocess
import os
import sys
import numpy as np
import update_settings as us
import parse_spe_reaction_info as psri
import prepare_path_name_time as ppnt
import parse_pattern as pp
import pattern_statistics as ps
import naming
import trajectory
import make_figures as mf
import global_settings


def quick_clean_up(data_dir, flag="", species_path=False):
    prefix = ""
    if species_path is True:
        prefix = "species_"

    if flag == "":
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate.csv")
    else:
        f_n_pn = os.path.join(data_dir, "output",
                              prefix + "pathway_name_candidate_" + str(flag) + ".csv")

    try:
        os.remove(f_n_pn)
    except OSError:
        pass

    if flag == "":
        f_n_pp = os.path.join(data_dir, "output",
                              prefix + "pathway_prob.csv")
    else:
        f_n_pp = os.path.join(data_dir, "output",
                              prefix + "pathway_prob_" + str(flag) + ".csv")

    try:
        os.remove(f_n_pp)
    except OSError:
        pass

    return


def update_terminal_species_setting(data_dir, terminal_spe=None):
    """
    update settings.json, primarily for terminal species
    """
    us.update_terminal_species_setting(data_dir, terminal_spe=terminal_spe)


def update_chattering_species_setting(data_dir, atom_followed="C"):
    """
    update settings.json, primarily for chattering species and fast reactions
    """
    us.update_chattering_species_setting(data_dir, atom_followed)


def copy_sohr_files(data_dir, species_path=False):
    """
    copy SOHR files from C++ routine
    """
    naming.copy_sohr_files(data_dir, species_path=species_path)


def species_count(data_dir, top_n=50, norm=False):
    """
    count species occurence in pathway, multiply by accurate pathwap probability
    """
    ps.species_count(data_dir, top_n=top_n, norm=norm)


def reaction_count(data_dir, top_n=50, norm=False):
    """
    count reaction occurence in pathway, multiply by accurate pathwap probability
    """
    ps.reaction_count(data_dir, top_n=top_n, norm=norm)


def initiation_reaction_count(data_dir, top_n=50, norm=False):
    """
    count initiation reaction occurence in pathway, multiply by accurate pathwap probability
    """
    ps.initiation_reaction_count(data_dir, top_n=top_n, norm=norm)


def species_cycle(data_dir, top_n=50, norm=False):
    """
    count species cycle in pathway, multiply by accurate pathwap probability
    """
    ps.species_cycle(data_dir, top_n=top_n, norm=norm)


def species_production_path(data_dir, spe='OH', top_n=50, norm=False):
    """
    count species production pathway or sub-pathway in pathway,
    multiply by accurate pathwap probability
    """
    ps.species_production_path(data_dir, spe=spe, top_n=top_n, norm=norm)


def species_production_reaction(data_dir, spe='OH', top_n=50, norm=False):
    """
    count species production reactions in pathway,
    multiply by accurate pathwap probability
    """
    ps.species_production_reaction(data_dir, spe=spe, top_n=top_n, norm=norm)


def symbolic_path_2_real_path(data_dir, top_n=50, flag="", end_s_idx=None, species_path=False, max_rows=5000):
    """
    convert symbolic pathway to real pathway with real species name and real reaction name
    flag indicates a specific job, for example, pathway end time = 1.0, the j-th run,
    any unique symbol shall work
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"

    if flag == "":
        out_file_name = prefix + "pathname_prob.csv"
    else:
        out_file_name = prefix + "pathname_prob_" + str(flag) + ".csv"

    path_stat_fn = prefix + "pathway_stat.csv"

    psri.symbolic_path_2_real_path(
        data_dir,
        os.path.join(
            data_dir, "output", path_stat_fn),
        os.path.join(
            data_dir, "output", out_file_name),
        top_n, end_s_idx, max_rows=max_rows)


def delete_non_dlsode_files(data_dir):
    """
    delete none dlsode files
    """
    os.chdir(data_dir)
    cmd = ["find", "./output", "-type", "f",
           "!", "-name", "*dlsode*", "-delete"]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

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


def make_run(src_dir, data_dir):
    """
    src_dir, source file (C++ exetutable file)
    make run
    """
    os.chdir(src_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

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


def run_dlsode(src_dir, data_dir, max_time, critical_time):
    """
    Run dlsode
    """
    os.chdir(src_dir)
    us.update_dlsode_setting(data_dir, max_time, critical_time)
    make_run(src_dir, data_dir)


def spe_concentration_at_time_w2f(src_dir, data_dir, tau=10.0, end_t=1.0):
    """
    write species concentration at a time to file
    """
    os.chdir(src_dir)
    us.update_spe_concentration_at_time_w2f(data_dir, tau=tau, end_t=end_t)
    make_run(src_dir, data_dir)


def run_mc_trajectory(src_dir, data_dir, n_traj=1000000, atom_followed="C", init_spe=114,
                      tau=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    Run mc trajectory
    """
    os.chdir(src_dir)
    us.update_mc_trajectory_setting(
        data_dir, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe,
        tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)
    make_run(src_dir, data_dir)


def evaluate_pathway_probability(
        src_dir, data_dir, top_n=5, num_t=1, flag="", n_traj=10000,
        atom_followed="C", init_spe=114, traj_max_t=100.0,
        tau=10.0, begin_t=0.0, end_t=1.0, top_n_s=10,
        spe_oriented=True, end_s_idx=None, species_path=False,
        path_reg=None, no_path_reg=None,
        spe_idx=None, spe_production_oriented=False,
        fixed_t0_or_tf=None, same_path_list=False):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(src_dir)

    if spe_oriented is True:
        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                data_dir, exclude=None, top_n=top_n_s, traj_max_t=traj_max_t,
                tau=tau, end_t=end_t, tag="M", atoms=[atom_followed])

        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path,
            path_reg=path_reg, no_path_reg=no_path_reg, spe_idx=spe_idx, spe_production_oriented=spe_production_oriented,
            same_path_list=same_path_list)

        ppnt.prepare_pathway_time(
            data_dir, top_n=n_path, num=num_t, flag=flag,
            begin_t=begin_t, end_t=end_t, species_path=species_path,
            fixed_t0_or_tf=fixed_t0_or_tf)

        us.update_eval_path_integral(
            data_dir, top_n=n_path, n_traj=n_traj,
            atom_followed=atom_followed, init_spe=init_spe,
            tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)

    else:
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path,
            path_reg=path_reg, no_path_reg=no_path_reg, spe_idx=spe_idx, spe_production_oriented=spe_production_oriented,
            same_path_list=same_path_list)
        ppnt.prepare_pathway_time(
            data_dir, top_n=n_path, num=num_t, flag=flag, begin_t=begin_t, end_t=end_t,
            species_path=species_path, fixed_t0_or_tf=fixed_t0_or_tf)

        us.update_eval_path_integral(
            data_dir, top_n=n_path, n_traj=n_traj, atom_followed=atom_followed, init_spe=init_spe,
            tau=tau, begin_t=begin_t, end_t=end_t, species_path=species_path)
    make_run(src_dir, data_dir)


def Merchant_f_2d_t0_tf(
        src_dir, data_dir, top_n=5, num_t=25, flag="", n_traj=10000,
        atom_followed="C", init_spe=114, traj_max_t=100.0,
        tau=10.0, begin_t=0.0, end_t=1.0,
        path_reg=None, no_path_reg=None,
        spe_idx=10, min_delta_t=None, num_delta_t=None, delta_t_vec=None):

    if flag == "":
        f_n_merchant_f = os.path.join(data_dir, "output",
                                      "Merchant_f_2d.csv")
    else:
        f_n_merchant_f = os.path.join(data_dir, "output",
                                      "Merchant_f_2d_" + str(flag) + ".csv")
    try:
        os.remove(f_n_merchant_f)
    except OSError:
        pass

    time_vec = np.linspace(begin_t, end_t, num_t)
    for i in range(num_t - 1):
        b_t = time_vec[i]
        if delta_t_vec is not None:
            end_t_vec = []
            for _, dt_val in enumerate(delta_t_vec):
                if b_t + float(dt_val) <= end_t:
                    end_t_vec.append(b_t + float(dt_val))
        elif min_delta_t is None or num_delta_t is None:
            end_t_vec = time_vec[i + 1:]
        else:
            end_t_vec = []
            for idx in range(int(num_delta_t)):
                if b_t + (idx + 1) * float(min_delta_t) <= end_t:
                    end_t_vec.append(b_t + (idx + 1) * float(min_delta_t))
        print(end_t_vec)
        for e_t in end_t_vec:
            # num_t set to be 1
            evaluate_pathway_probability(
                src_dir, data_dir, top_n=top_n, num_t=1, flag=flag, n_traj=n_traj,
                atom_followed=atom_followed, init_spe=init_spe, traj_max_t=traj_max_t,
                tau=tau, begin_t=b_t, end_t=e_t, top_n_s=None,
                spe_oriented=False, end_s_idx=None, species_path=False,
                path_reg=path_reg, no_path_reg=no_path_reg,
                spe_idx=spe_idx, spe_production_oriented=True,
                fixed_t0_or_tf="t0", same_path_list=True)

            f_n_pp = os.path.join(data_dir, "output",
                                  "pathway_prob" + ".csv")
            f_vec = np.loadtxt(f_n_pp, dtype=float, delimiter=',')

            text = str(b_t) + ',' + str(e_t)
            for _, val in enumerate(f_vec):
                text += (',' + str(val))

            f_sum = np.sum(f_vec)
            text += (',' + str(f_sum) + '\n')
            with open(f_n_merchant_f, 'a') as f_handler:
                f_handler.write(text)

    return


def spe_concentration_converge_at_different_times(
        src_dir, data_dir, flag="", top_n=5,
        mc_n_traj=10000, pi_n_traj=10000,
        traj_max_t=100.0, tau=10.0,
        atom_followed="C",
        init_spe=114,
        path_reg=None, no_path_reg=None,
        species_path=False,
        begin_t=0.0,
        end_t_vec=[0.1, 0.2]):
    """
    For a single species, at different time points, run monte carlo pathway generation,
    select top n pathways, evaluate pathway probabilities
    species concentration converge --> scc
    end species index control by regular expression "path_reg"
    """
    prefix = ""
    if species_path is True:
        prefix = "species_"

    # save end time vector
    if flag == "":
        f_n_scc_time = os.path.join(data_dir, "output",
                                    prefix + "scc_time.csv")
    else:
        f_n_scc_time = os.path.join(data_dir, "output",
                                    prefix + "scc_time_" + str(flag) + ".csv")
    try:
        os.remove(f_n_scc_time)
    except OSError:
        pass
    np.savetxt(f_n_scc_time, end_t_vec, fmt='%.18e', delimiter='\n')

    # pathway probabilities
    if flag == "":
        f_n_scc_pp = os.path.join(data_dir, "output",
                                  prefix + "scc_path_prob.csv")
    else:
        f_n_scc_pp = os.path.join(data_dir, "output",
                                  prefix + "scc_path_prob_" + str(flag) + ".csv")
    try:
        os.remove(f_n_scc_pp)
    except OSError:
        pass
    # temporary files
    f_n_pp = os.path.join(data_dir, "output",
                          prefix + "pathway_prob" + ".csv")

    for e_t in end_t_vec:
        print("time point:\t", e_t)
        # run mc trajectory
        run_mc_trajectory(src_dir, data_dir, n_traj=mc_n_traj, atom_followed=atom_followed, init_spe=init_spe,
                          tau=tau, begin_t=begin_t, end_t=e_t, species_path=species_path)
        # prepare pathway name candidate and evaluate pathway probabilities
        evaluate_pathway_probability(
            src_dir, data_dir, top_n=top_n, num_t=1, flag=flag, n_traj=pi_n_traj,
            atom_followed=atom_followed, init_spe=init_spe, traj_max_t=traj_max_t,
            tau=tau, begin_t=begin_t, end_t=e_t, top_n_s=None,
            spe_oriented=False, end_s_idx=None, species_path=species_path,
            path_reg=path_reg, no_path_reg=no_path_reg,
            spe_idx=None, spe_production_oriented=False,
            fixed_t0_or_tf="t0", same_path_list=False)
        # save pathway probabilities at current time to a file
        data_pp = np.loadtxt(f_n_pp, dtype=float, delimiter=",")
        with open(f_n_scc_pp, 'a') as f_handler:
            np.savetxt(f_handler, data_pp, fmt='%.18e',
                       delimiter=',', newline='\n', header='')
        return


def evaluate_pathway_AT(src_dir, data_dir, top_n=5, flag="", n_traj=10000,
                        traj_max_t=100.0, atom_followed="C",
                        tau=10.0, begin_t=0.0, end_t=1.0,
                        top_n_s=10, spe_oriented=True, end_s_idx=None, species_path=False):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(src_dir)

    if spe_oriented is True:
        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                data_dir, exclude=None, top_n=top_n_s, traj_max_t=traj_max_t,
                tau=tau, end_t=end_t, tag="M", atoms=[atom_followed])
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag,
            end_s_idx=end_s_idx,
            species_path=species_path,
            same_path_list=True)

        us.update_eval_path_AT(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)

    else:
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag,
            end_s_idx=None,
            species_path=species_path,
            same_path_list=True)
        us.update_eval_path_AT(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)

    make_run(src_dir, data_dir)


def evaluate_pathway_AT_no_IT(src_dir, data_dir, top_n=5, flag="", n_traj=10000,
                              atom_followed="C", traj_max_t=100.0,
                              tau=10.0, begin_t=0.0, end_t=1.0, top_n_s=10,
                              spe_oriented=True, end_s_idx=None, species_path=False):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(src_dir)

    if spe_oriented is True:
        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                data_dir, exclude=None, top_n=top_n_s, traj_max_t=traj_max_t,
                tau=tau, end_t=end_t, tag="M", atoms=[atom_followed])
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path)

        us.update_eval_path_AT_no_IT(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)
    else:
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=None, species_path=species_path)
        us.update_eval_path_AT_no_IT(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)

    make_run(src_dir, data_dir)


def evaluate_pathway_AT_with_SP(src_dir, data_dir, top_n=5, flag="", n_traj=10000,
                                atom_followed="C", traj_max_t=100.0,
                                tau=10.0, begin_t=0.0, end_t=1.0,
                                top_n_s=10, spe_oriented=True, end_s_idx=None, species_path=False):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(src_dir)

    if spe_oriented is True:
        if end_s_idx is None or end_s_idx is []:
            end_s_idx, _, _ = trajectory.get_species_with_top_n_concentration(
                data_dir, exclude=None, top_n=top_n_s, traj_max_t=traj_max_t,
                tau=tau, end_t=end_t, tag="M", atoms=[atom_followed])
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=end_s_idx, species_path=species_path)

        us.update_eval_path_AT_with_SP(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)
    else:
        n_path = ppnt.prepare_pathway_name(
            data_dir, top_n=top_n, flag=flag, end_s_idx=None, species_path=species_path)
        us.update_eval_path_AT_with_SP(
            data_dir, top_n=n_path, n_traj=n_traj,
            tau=tau, begin_t=begin_t, end_t=end_t)

    make_run(src_dir, data_dir)


def evaluate_passage_time_of_species(src_dir, data_dir, flag="", n_traj=10000,
                                     tau=10.0, begin_t=0.0, end_t=1.0, init_s_idx=None):
    """
    evaluate pathway probability
    top_n_s is top N species number
    num_t is number of time points
    """
    os.chdir(data_dir)

    # species_path have to be true, just how the C++ codes written
    # physically, passage time don't depend on reaction details
    # only depends on species survial probability
    n_path = ppnt.prepare_pathway_name_for_passage_time(
        data_dir, flag=flag, init_s_idx=init_s_idx)
    us.update_eval_path_AT(
        data_dir, top_n=n_path, n_traj=n_traj,
        tau=tau, begin_t=begin_t, end_t=end_t)

    make_run(src_dir, data_dir)


# http://stackoverflow.com/questions/3000724/running-matlab-in-the-background
def make_a_figure(data_dir, ind):
    """
    make a figure
    """
    os.chdir(data_dir)
    matlab_script_dir = os.path.join(
        data_dir, "tools/data_analysis/H2_O2_reaction_network_ODE_solver")
    matlab_script_filename = "plot_concentration_v2"
    matlab_cmd = "cd " + matlab_script_dir + "; " + matlab_script_filename + "(" + str(
        ind) + ")" + "; cd ../../..; exit;"
    cmd = ["nohup", "matlab", "-nosplash", "-nodisplay", "-r", matlab_cmd]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

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


def make_figures(data_dir):
    """
    make figures
    """
    for i in range(1, 9):
        make_a_figure(data_dir, i)


def propane_make_figures(data_dir, species_path=False):
    """
    make figures for propane system
    """
    g_s = global_settings.get_setting(data_dir)

    spe_idx, _, spe_exclude_name = trajectory.get_species_with_top_n_concentration(
        data_dir, exclude=None, top_n=g_s['top_n_s'], traj_max_t=g_s['traj_max_t'],
        tau=g_s['tau'], end_t=g_s['end_t'], tag=g_s['tag'], atoms=[g_s['atom_f']])
    mf.plot_concentrations(
        data_dir, spe_idx=spe_idx, tau=g_s['tau'], end_t=g_s['end_t'], tag=g_s['tag'],
        exclude_names=spe_exclude_name, renormalization=True)
    mf.plot_reaction_rates(
        data_dir, reaction_idx=[1068, 1070, 1072, 1074, 1076], tau=g_s['tau'], end_t=1.0, tag=g_s['tag'])
    for s_i in spe_idx:
        print(spe_idx)
        mf.plot_species_pathway_prob(data_dir, top_n=g_s['top_n_p'],
                                     exclude_names=spe_exclude_name,
                                     init_spe=g_s['init_s'], atom_followed=g_s['atom_f'],
                                     tau=g_s['tau'], end_t=g_s['end_t'],
                                     end_s_idx=s_i, species_path=species_path)
    mf.plot_reaction_rate_constant(data_dir)
    mf.plot_top_n_spe_concentration(
        data_dir, exclude_names=None, atom_followed=g_s['atom_f'], end_t=g_s['end_t'], top_n=10)


def send_email(data_dir):
    """
    send email to elliot.srbai@gmail.com
    """
    os.chdir(data_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
           "-u", "RUNNING JOB", "-m", "JOB FINISHED." + "\n" + data_dir,
           "-a", os.path.join(data_dir, "output", "output_all.txt")]

    # Open/Create the output file
    out_file = open(os.path.join(
        data_dir, 'output', 'output_all.txt'), 'ab+')
    error_file = open(os.path.join(
        data_dir, 'output', 'error_all.txt'), 'ab+')

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


if __name__ == '__main__':
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print("test")
    # symbolic_path_2_real_path_pff(DATA_DIR, "heuristic_pathname_O_10_10_3.csv")
