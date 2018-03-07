"""
main driver
"""

import os
import sys
import time

import update_settings as us
import job_drivers
import global_settings
import pattern_statistics as ps

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
    # job_drivers.run_dlsode(DATA_DIR, G_S['traj_max_t'], G_S['traj_critical_t'])

    # update terminal species
    job_drivers.update_terminal_species_setting(DATA_DIR, G_S['terminal_spe'])

    # update chattering species and fast reactions
    job_drivers.update_chattering_species_setting(
        DATA_DIR, G_S['atom_f'])

    # run monte carlo trajectory
    job_drivers.run_mc_trajectory(
        SRC_DIR, DATA_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
        init_spe=G_S['init_s'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['mc_t'],
        species_path=G_S['species_path'])

    # evaluate path integral-->pathway probability
    job_drivers.evaluate_pathway_probability(
        SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
        n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
        traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
        top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
        end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # job_drivers.evaluate_pathway_AT(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # job_drivers.evaluate_pathway_AT_no_IT(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # job_drivers.evaluate_pathway_AT_with_SP(
    #     SRC_DIR, DATA_DIR, top_n=G_S['top_n_p'], flag="",
    #     n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
    #     traj_max_t=G_S['traj_max_t'], tau=G_S['tau'], begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #     top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'],
    #     end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # job_drivers.evaluate_passage_time_of_species(SRC_DIR, DATA_DIR, flag="", n_traj=G_S['pi_n_traj'],
    #                                              atom_followed=G_S['atom_f'], init_spe=G_S['init_s'], tau=G_S['tau'],
    #                                              begin_t=G_S['begin_t'], end_t=G_S['end_t'],
    #                                              end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # convert symbolic pathway to real pathway
    # with real species names and real reaction expression
    job_drivers.symbolic_path_2_real_path(DATA_DIR, top_n=G_S['top_n_p'], flag="",
                                          end_s_idx=None, species_path=G_S['species_path'])

    # copy SOHR/C++ routine files
    job_drivers.copy_sohr_files(DATA_DIR, species_path=G_S['species_path'])

    # ps.parse_spe_production_along_path(DATA_DIR, top_n=G_S['top_n_p'], spe_idx=10,
    #                                    init_spe=G_S['init_s'], atom_followed=G_S['atom_f'],
    #                                    end_t=G_S['end_t'], species_path=G_S['species_path'],
    #                                    axis=0, path_branching_factor=False)

    # # # species count
    # # job_drivers.species_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # reaction count
    # # job_drivers.reaction_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # initiation reaction count
    # # job_drivers.initiation_reaction_count(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species cycle
    # # job_drivers.species_cycle(DATA_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species production path
    # # job_drivers.species_production_path(
    # #     DATA_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # # # species production reaction
    # # job_drivers.species_production_reaction(
    # #     DATA_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # # propane make figures
    # job_drivers.propane_make_figures(
    #    DATA_DIR, species_path=G_S['species_path'])

    # send email
    # job_drivers.send_email(DATA_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
