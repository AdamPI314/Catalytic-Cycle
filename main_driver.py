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

    # run monte carlo trajectory
    j_b.run_mc_trajectory(
        SRC_DIR, DATA_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
        init_spe=G_S['init_s'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['mc_t'],
        species_path=G_S['species_path'])

    # # evaluate path integral-->pathway probability
    # j_b.evaluate_pathway_probability(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'], path_reg='^S62R(736|738)')

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
    #     SRC_DIR, DATA_DIR, flag="", n_traj=G_S['pi_n_traj'],
    #     tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     init_s_idx=G_S['end_s_idx'])

    # # convert symbolic pathway to real pathway
    # # with real species names and real reaction expression
    # j_b.symbolic_path_2_real_path(DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #                               end_s_idx=None, species_path=G_S['species_path'])

    # if G_S['species_path'] is False:
    #     psri.symbolic_path_2_real_path_pff(
    #         DATA_DIR, 'pathway_name_candidate.csv')
    # else:
    #     psri.symbolic_path_2_real_path_pff(
    #         DATA_DIR, 'species_pathway_name_candidate.csv')

    # copy SOHR/C++ routine files
    j_b.copy_sohr_files(DATA_DIR, species_path=G_S['species_path'])

    # ps.parse_spe_production_along_path(
    #     DATA_DIR, top_n=G_S['top_n_p'], spe_idx=10,
    #     init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #     end_t=G_S['end_t'], species_path=G_S['species_path'],
    #     axis=0, path_branching_factor=False)

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
