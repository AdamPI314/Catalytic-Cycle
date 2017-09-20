"""
update settings.json
"""

import os
from shutil import copy2
import read_write_configuration as rwc
import global_settings


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
    setting['time']['tau'] = time

    setting['job']['job_type'] = "solve_ODEs_for_concentration_using_LSODE"
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_spe_concentration_at_time_w2f(file_dir, tau=1.0):
    """
    update settings.json, primarily for update_spe_concentration_at_time_w2f
    """
    # there will always be a current setting
    fn0 = os.path.join(file_dir, "input", "setting_backup.json")
    fn1 = os.path.join(file_dir, "input", "setting.json")

    if not os.path.isfile(fn0):
        copy2(fn1, fn0)

    setting = rwc.read_configuration(
        os.path.join(file_dir, 'input', 'setting.json'))

    setting['job']['job_type'] = "write_concentration_at_time_to_file"
    setting['pathway']['tau'] = tau

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))


def update_mc_trajectory_setting(file_dir, time=1.0, n_traj=1000000, atom_followed="C", init_spe=114, max_tau=1.0):
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

    fast_reaction, trapped_spe = global_settings.get_fast_rxn_trapped_spe(file_dir)
    setting['pathway']['fast_reaction'] = fast_reaction
    setting['pathway']['trapped_species'] = trapped_spe

    setting['time']['critical_time'] = time
    setting['time']['max_time'] = time
    setting['time']['tau'] = time
    setting['pathway']['trajectoryNumber'] = n_traj
    setting['pathway']['atom_followed'] = atom_followed
    setting['pathway']['init_spe'] = init_spe
    setting['pathway']['tau'] = max_tau

    setting['job']['job_type'] = "generate_pathway_running_Monte_carlo_trajectory"
    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
    return


def update_eval_path_integral(file_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114, max_tau=1.0):
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

    fast_reaction, trapped_spe = global_settings.get_fast_rxn_trapped_spe(file_dir)
    setting['pathway']['fast_reaction'] = fast_reaction
    setting['pathway']['trapped_species'] = trapped_spe

    setting['job']['job_type'] = "evaluate_path_integral_over_time"
    setting['pathway']['pathwayEndWith'] = "ALL"
    setting['pathway']['topN'] = [top_n]
    setting['pathway']['trajectoryNumber'] = n_traj
    setting['pathway']['atom_followed'] = atom_followed
    setting['pathway']['init_spe'] = init_spe
    setting['pathway']['tau'] = max_tau

    rwc.write_configuration(setting, os.path.join(
        file_dir, 'input', 'setting.json'))
