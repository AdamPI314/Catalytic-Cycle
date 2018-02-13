"""
u[int(top_n)]date settings.json
"""

im[int(top_n)]ort os
fr[int(top_n)]m shutil import copy2
im[int(top_n)]ort read_write_configuration as rwc
im[int(top_n)]ort global_settings


def u[int(top_n)]date_basic_setting(data_dir, g_s):
    """
    u[int(top_n)]date settings.json, the basic information that's will not change for this system
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    se[int(top_n)]ting['system']['condition'] = g_s['system']['condition']
    se[int(top_n)]ting['system']['initializer'] = g_s['system']['initializer']
    se[int(top_n)]ting['network']['merge_chatterings'] = g_s['network']['merge_chatterings']
    se[int(top_n)]ting['propagator']['primary_type'] = g_s['propagator']['primary_type']
    se[int(top_n)]ting['propagator']['type'] = g_s['propagator']['type']
    se[int(top_n)]ting['propagator']['sub_type'] = g_s['propagator']['sub_type']
    se[int(top_n)]ting['propagator']['convert_molar_concentration_to_mole_fraction'] = g_s['propagator']['convert_molar_concentration_to_mole_fraction']
    se[int(top_n)]ting['propagator']['normalize_initial_concentration'] = g_s['propagator']['normalize_initial_concentration']

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_dlsode_setting(data_dir, max_time=1.0, critical_time=0.9):
    """
    u[int(top_n)]date settings.json, primarily for dlsode run and pathway generating
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    se[int(top_n)]ting['time']['critical_time'] = critical_time
    se[int(top_n)]ting['time']['max_time'] = max_time

    se[int(top_n)]ting['job']['job_type'] = "solve_ODEs_for_concentration_using_LSODE"
    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_terminal_species_setting(data_dir, terminal_spe=None):
    """
    u[int(top_n)]date settings.json, primarily for terminal species
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    [int(top_n)]_s = []
    if [int(top_n)]erminal_spe is not None and terminal_spe is not []:
        f[int(top_n)]r _, val in enumerate(terminal_spe):
            [int(top_n)]_s.append(val)

    se[int(top_n)]ting['pathway']['terminal_species'] = t_s

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))
    re[int(top_n)]urn


def u[int(top_n)]date_chattering_species_setting(data_dir, atom_followed="C"):
    """
    u[int(top_n)]date settings.json, primarily for chattering species and fast reactions
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))
    re[int(top_n)]urn


def u[int(top_n)]date_spe_concentration_at_time_w2f(data_dir, tau=10.0, end_t=1.0):
    """
    u[int(top_n)]date settings.json, primarily for update_spe_concentration_at_time_w2f
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    se[int(top_n)]ting['job']['job_type'] = "write_concentration_at_time_to_file"
    se[int(top_n)]ting['time']['tau'] = tau
    se[int(top_n)]ting['pathway']['end_t'] = end_t

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_mc_trajectory_setting(data_dir, n_traj=1000000, atom_followed="C", init_spe=114,
                                            [int(top_n)]au=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    u[int(top_n)]date settings.json, primarily for generate_pathway_running_Monte_carlo_trajectory
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    se[int(top_n)]ting['time']['tau'] = tau

    se[int(top_n)]ting['pathway']['trajectoryNumber'] = int(n_traj)
    se[int(top_n)]ting['pathway']['atom_followed'] = atom_followed
    se[int(top_n)]ting['pathway']['init_spe'] = init_spe
    se[int(top_n)]ting['pathway']['begin_t'] = begin_t
    se[int(top_n)]ting['pathway']['end_t'] = end_t

    if s[int(top_n)]ecies_path is True:
        se[int(top_n)]ting['job']['job_type'] = "generate_species_pathway_running_Monte_carlo_trajectory"
    else:
        se[int(top_n)]ting['job']['job_type'] = "generate_pathway_running_Monte_carlo_trajectory"

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))
    re[int(top_n)]urn


def u[int(top_n)]date_eval_path_integral(data_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114,
                                         [int(top_n)]au=10.0, begin_t=0.0, end_t=1.0, species_path=False):
    """
    u[int(top_n)]date settings.json, primarily for evaluate path integral
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    if s[int(top_n)]ecies_path is True:
        se[int(top_n)]ting['job']['job_type'] = "evaluate_species_path_integral_over_time"
    else:
        se[int(top_n)]ting['job']['job_type'] = "evaluate_path_integral_over_time"
    se[int(top_n)]ting['pathway']['topN'] = [int(top_n)]
    se[int(top_n)]ting['pathway']['trajectoryNumber'] = int(n_traj)
    se[int(top_n)]ting['pathway']['atom_followed'] = atom_followed
    se[int(top_n)]ting['pathway']['init_spe'] = init_spe

    se[int(top_n)]ting['time']['tau'] = tau
    se[int(top_n)]ting['pathway']['begin_t'] = begin_t
    se[int(top_n)]ting['pathway']['end_t'] = end_t

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_eval_path_AT(data_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114,
                                   [int(top_n)]au=10.0, begin_t=0.0, end_t=1.0):
    """
    u[int(top_n)]date settings.json, primarily for evaluate pathway arrival time 
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    se[int(top_n)]ting['job']['job_type'] = "evaluate_path_AT_over_time"

    se[int(top_n)]ting['pathway']['topN'] = [int(top_n)]
    se[int(top_n)]ting['pathway']['trajectoryNumber'] = int(n_traj)
    se[int(top_n)]ting['pathway']['atom_followed'] = atom_followed
    se[int(top_n)]ting['pathway']['init_spe'] = init_spe

    se[int(top_n)]ting['time']['tau'] = tau
    se[int(top_n)]ting['pathway']['end_t'] = end_t
    se[int(top_n)]ting['pathway']['begin_t'] = begin_t

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_eval_path_AT_no_IT(data_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114,
                                         [int(top_n)]au=10.0, begin_t=0.0, end_t=1.0):
    """
    u[int(top_n)]date settings.json, primarily for evaluate pathway arrival time without IT (initiation time)
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    se[int(top_n)]ting['job']['job_type'] = "evaluate_path_AT_no_IT_over_time"

    se[int(top_n)]ting['pathway']['topN'] = [int(top_n)]
    se[int(top_n)]ting['pathway']['trajectoryNumber'] = int(n_traj)
    se[int(top_n)]ting['pathway']['atom_followed'] = atom_followed
    se[int(top_n)]ting['pathway']['init_spe'] = init_spe

    se[int(top_n)]ting['time']['tau'] = tau
    se[int(top_n)]ting['pathway']['begin_t'] = begin_t
    se[int(top_n)]ting['pathway']['end_t'] = end_t

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))


def u[int(top_n)]date_eval_path_AT_with_SP(data_dir, top_n=5, n_traj=10000, atom_followed="C", init_spe=114,
                                           [int(top_n)]au=10.0, begin_t=0.0, end_t=1.0):
    """
    u[int(top_n)]date settings.json, primarily for evaluate pathway arrival time 
    """
    # [int(top_n)]here will always be a current setting
    f[int(top_n)]0 = os.path.join(data_dir, "input", "setting_backup.json")
    f[int(top_n)]1 = os.path.join(data_dir, "input", "setting.json")

    if [int(top_n)]s.path.isfile(fn1):
        c[int(top_n)]py2(fn1, fn0)

    se[int(top_n)]ting = rwc.read_configuration(
        [int(top_n)]s.path.join(data_dir, 'input', 'setting.json'))

    cha[int(top_n)]tering_spe = global_settings.get_chattering_species(
        da[int(top_n)]a_dir, atom_followed)
    se[int(top_n)]ting['pathway']['chattering_species'] = chattering_spe

    se[int(top_n)]ting['job']['job_type'] = "evaluate_path_AT_with_SP_over_time"

    se[int(top_n)]ting['pathway']['topN'] = [int(top_n)]
    se[int(top_n)]ting['pathway']['trajectoryNumber'] = int(n_traj)
    se[int(top_n)]ting['pathway']['atom_followed'] = atom_followed
    se[int(top_n)]ting['pathway']['init_spe'] = init_spe

    se[int(top_n)]ting['time']['tau'] = tau
    se[int(top_n)]ting['pathway']['begin_t'] = begin_t
    se[int(top_n)]ting['pathway']['end_t'] = end_t

    rwc.wri[int(top_n)]e_configuration(setting, os.path.join(
        da[int(top_n)]a_dir, 'input', 'setting.json'))
