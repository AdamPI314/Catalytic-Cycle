"""
update settings.json
"""

import os
from shutil import copy2
import read_write_configuration as rwc


def update_dlsode_setting(file_dir, time=1.0):
    """
    update settings.json, primarily for dlsode run and pathway generating
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if not os.path.isfile(fn0):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['time']['critical_time'] = time
    setting['time']['max_time'] = time
    setting['time']['path_end_time'] = time

    setting['job']['job_type'] = "solve_ODEs_for_concentration_using_LSODE"
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_mc_trajectory_setting(file_dir, time=1.0):
    """
    update settings.json, primarily for generate_pathway_running_Monte_carlo_trajectory
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if not os.path.isfile(fn0):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['time']['critical_time'] = time
    setting['time']['max_time'] = time
    setting['time']['path_end_time'] = time
    setting['pathway']['trajectoryNumber'] = 1000000

    setting['job']['job_type'] = "generate_pathway_running_Monte_carlo_trajectory"
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_eval_path_integral(file_dir, top_n=5):
    """
    update settings.json, primarily for evaluate path integral
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if not os.path.isfile(fn0):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['job']['job_type'] = "evaluate_path_integral_over_time"
    setting['pathway']['pathwayEndWith'] = "ALL"
    setting['pathway']['topN'] = [top_n]
    setting['pathway']['trajectoryNumber'] = 10000
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
