"""
Job drivers
"""

import subprocess
import os
import update_settings as us
import parse_spe_reaction_info as psri


def symbolic_path_2_real_path(file_dir, top_n=50, flag=""):
    """
    convert symbolic pathway to real pathway with real species name and real reactin name
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
        top_n)


def delete_non_dlsode_files(file_dir):
    """
    delete none dlsode files
    """
    os.chdir(file_dir)
    cmd = ["find", "./output", "-type", "f",
           "!", "-name", "*dlsode*", "-delete"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'a+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'a+')

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


def path_ode_run(file_dir):
    """
    make run
    """
    os.chdir(file_dir)
    cmd = ["make", "run"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'a+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'a+')

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


def run_dlsode(file_dir):
    """
    Run dlsode
    """
    os.chdir(file_dir)
    us.update_dlsode_settings(file_dir, time=None)
    path_ode_run(file_dir)

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
        file_dir, 'output', 'output_all.txt'), 'a+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'a+')

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


def send_email(file_dir):
    """
    send email to elliot.srbai@gmail.com
    """
    os.chdir(file_dir)
    cmd = ["sendemail", "-f", "elliot.srbai@gmail.com", "-t", "bunnysirah@hotmail.com",
           "-u", "RUNNING JOB", "-m", "JOB FINISHED." + "\n" + file_dir,
           "-a", "./output/general_output.out"]

    # Open/Create the output file
    out_file = open(os.path.join(
        file_dir, 'output', 'output_all.txt'), 'a+')
    error_file = open(os.path.join(
        file_dir, 'output', 'error_all.txt'), 'a+')

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
