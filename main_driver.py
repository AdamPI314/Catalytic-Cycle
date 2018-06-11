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

    # run monte carlo trajectory
    j_b.run_mc_trajectory(
        SRC_DIR, DATA_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
        init_spe=G_S['init_s'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['mc_t'],
        species_path=G_S['species_path'])

    # evaluate path integral-->pathway probability, multiple use
    j_b.evaluate_pathway_probability(
        SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
        n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
        traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
        top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
        end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'],
        path_reg=G_S['path_reg'], no_path_reg=G_S['no_path_reg'],
        spe_idx=None, spe_production_oriented=True,
        fixed_t0_or_tf=G_S['fixed_t0_or_tf'])

    # j_b.Merchant_f_2d_t0_tf(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=25, flag="", n_traj=G_S['pi_n_traj'],
    #     atom_followed=G_S['atom_f'], init_spe=G_S['init_s'], traj_max_t=G_S['traj_max_t'],
    #     tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     path_reg=G_S['path_reg'], no_path_reg=G_S['no_path_reg'],
    #     spe_idx=10, min_delta_t=None, num_delta_t=None,
    #     delta_t_vec=[1.2859087486111409e-06, 3.214771871527852e-06, 6.429543743055704e-06, 9.644315614583557e-06, 1.285908748611141e-05, 3.2147718715278524e-05, 6.429543743055705e-05, 9.644315614583559e-05, 0.0001285908748611141, 0.00032147718715278527, 0.0006429543743055705, 0.0009644315614583557, 0.001285908748611141, 0.0032147718715278524, 0.006429543743055705, 0.009644315614583556, 0.01285908748611141, 0.03214771871527852, 0.06429543743055705, 0.09644315614583557, 0.1285908748611141, 0.3214771871527852, 0.6429543743055705, 0.9644315614583557])

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
    #     mc_end_t_threshlod=0.15,
    #     end_t_vec=[0.0482215780729, 0.0964431561458, 0.144664734219, 0.192886312292, 0.196270282683, 0.199654253074, 0.203038223465, 0.206422193856, 0.209806164247, 0.213190134638, 0.216574105029, 0.21995807542, 0.223342045811, 0.226726016202, 0.230109986594, 0.233493956985, 0.236877927376, 0.240261897767, 0.243645868158, 0.247029838549, 0.25041380894, 0.253797779331, 0.257181749722])

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

    # copy SOHR/C++ routine files
    j_b.copy_sohr_files(DATA_DIR, species_path=G_S['species_path'])

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
