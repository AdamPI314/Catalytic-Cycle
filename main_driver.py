"""
main driver
"""

import os
import sys
import time

import update_settings as us
import job_drivers as j_b
import global_settings
# import pattern_statistics as ps
import parse_spe_reaction_info as psri
# import trajectory as traj

if __name__ == '__main__':
    TIME_I = time.time()

    # source file directory
    SRC_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(SRC_DIR)
    # data directory
    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))
    print(DATA_DIR)

    G_S = global_settings.get_setting(DATA_DIR)
    us.update_basic_setting(DATA_DIR, G_S)

    # # run dlosde
    # j_b.run_dlsode(DATA_DIR, G_S['traj_max_t'], G_S['traj_critical_t'])

    # update terminal species
    j_b.update_terminal_species_setting(DATA_DIR, G_S['terminal_spe'])

    # update chattering species and fast reactions
    j_b.update_chattering_species_setting(
        DATA_DIR, G_S['atom_f'])

    # quick clean up, remove a few files
    j_b.quick_clean_up(DATA_DIR, flag="", species_path=G_S['species_path'])

    # # run monte carlo trajectory
    # j_b.run_mc_trajectory(
    #     SRC_DIR, DATA_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
    #     init_spe=G_S['init_s'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['mc_t'],
    #     species_path=G_S['species_path'])

    # # evaluate path integral-->pathway probability, multiple use
    # j_b.evaluate_pathway_probability(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'],
    #     path_reg=G_S['path_reg'], no_path_reg=G_S['no_path_reg'],
    #     spe_idx=None, spe_production_oriented=True,
    #     fixed_t0_or_tf=G_S['fixed_t0_or_tf'],
    #     same_path_list=True)

    j_b.Merchant_f_2d_t0_tf(
        SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=10, flag="",
        mc_n_traj=G_S['mc_n_traj'], n_traj=G_S['pi_n_traj'],
        atom_followed=G_S['atom_f'], init_spe=G_S['init_s'], traj_max_t=G_S['traj_max_t'],
        tau=G_S['tau'], mc_end_t=G_S['mc_t'], t0_min=0.0, t0_max=0.9, end_t=G_S['end_t'],
        path_reg=G_S['path_reg'], no_path_reg=G_S['no_path_reg'],
        spe_idx=10, min_delta_t=None, num_delta_t=None,
        delta_t_vec=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9])

    # j_b.spe_concentration_converge_at_different_times(
    #     SRC_DIR, DATA_DIR, flag="",
    #     top_n=G_S['top_n_p'],
    #     mc_n_traj=G_S['mc_n_traj'], pi_n_traj=G_S['pi_n_traj'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'],
    #     atom_followed=G_S['atom_f'],
    #     init_spe=G_S['init_s'],
    #     path_reg=G_S['path_reg'],
    #     no_path_reg=G_S['no_path_reg'],
    #     species_path=G_S['species_path'],
    #     begin_t=0.0,
    #     mc_end_t_threshlod=0.25,
    #     end_t_vec=[0.0107159062384, 0.0214318124768, 0.0321477187152, 0.0428636249537, 0.0535795311921, 0.0642954374305, 0.0750113436689, 0.0857272499073, 0.0964431561457, 0.107159062384, 0.117874968623, 0.128590874861, 0.139306781099, 0.150022687338, 0.160738593576, 0.171454499815, 0.182170406053, 0.192886312291, 0.20360221853, 0.214318124768, 0.225034031007, 0.235749937245, 0.246465843484, 0.257181749722])

    # j_b.evaluate_pathway_AT(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # j_b.evaluate_pathway_AT_no_IT(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # j_b.evaluate_pathway_AT_with_SP(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # j_b.evaluate_passage_time_of_species(
    #     SRC_DIR, DATA_DIR, flag="", n_traj=1000000,
    #     tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     init_s_idx=G_S['init_s_idx'])

    # traj.cal_passage_time_distribution(
    #     DATA_DIR, G_S['init_s_idx'][0], G_S['tau'], t_f=G_S['end_t'], n_point=100000)

    # # convert symbolic pathway to real pathway
    # # with real species names and real reaction expression
    # j_b.symbolic_path_2_real_path(DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #                               end_s_idx=None, species_path=G_S['species_path'])

    if G_S['species_path'] is False:
        psri.symbolic_path_2_real_path_pff(
            DATA_DIR, 'pathway_name_candidate.csv')
    else:
        psri.symbolic_path_2_real_path_pff(
            DATA_DIR, 'species_pathway_name_candidate.csv')

    # # copy SOHR/C++ routine files
    # j_b.copy_sohr_files(DATA_DIR, species_path=G_S['species_path'])

    # ps.parse_spe_production_along_path(
    #     DATA_DIR, top_n=G_S['top_n_p'], spe_idx=[10],
    #     init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #     end_t=G_S['end_t'], species_path=G_S['species_path'],
    #     axis=0, path_branching_factor=False,
    #     s_consumption=False, s_production=True)

    # # # species count
    # # j_b.species_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # reaction count
    # # j_b.reaction_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # initiation reaction count
    # # j_b.initiation_reaction_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species cycle
    # # j_b.species_cycle(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species production path
    # # j_b.species_production_path(
    # #     DATA_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # # # species production reaction
    # # j_b.species_production_reaction(
    # #     DATA_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # # propane make figures
    # j_b.propane_make_figures(
    #    DATA_DIR, species_path=G_S['species_path'])

    # send email
    # j_b.send_email(DATA_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
