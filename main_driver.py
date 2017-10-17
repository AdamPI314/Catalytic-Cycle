"""
main driver
"""

import os
import sys
import time

import job_drivers
import global_settings

if __name__ == '__main__':
    TIME_I = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    G_S = global_settings.get_setting()

    # # run dlosde
    # job_drivers.run_dlsode(FILE_DIR, G_S['end_t'])

    # update trapped species and fast reactions
    job_drivers.update_trapped_species_fast_reaction_setting(FILE_DIR)

    # write specie concentration at a time to file
    job_drivers.spe_concentration_at_time_w2f(
        FILE_DIR, max_tau=G_S['max_tau'], tau=G_S['tau'])

    # run monte carlo trajectory
    job_drivers.run_mc_trajectory(
        FILE_DIR, n_traj=G_S['mc_n_traj'], atom_followed=G_S['atom_f'],
        init_spe=G_S['init_s'], max_tau=G_S['max_tau'], tau=G_S['tau'], species_path=G_S['species_path'])

    # evaluate path integral-->pathway probability
    job_drivers.evaluate_pathway_probability(
        FILE_DIR, top_n=G_S['top_n_p'], num_t=G_S['pi_n_time'], flag="",
        n_traj=G_S['pi_n_traj'], atom_followed=G_S['atom_f'], init_spe=G_S['init_s'],
        traj_end_time=G_S['end_t'], max_tau=G_S['max_tau'], tau=G_S['tau'],
        top_n_s=G_S['top_n_s'], spe_oriented=G_S['spe_oriented'], end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])

    # convert symbolic pathway to real pathway
    # with real species names and real reaction expression
    if G_S['spe_oriented'] is True:
        job_drivers.symbolic_path_2_real_path(
            FILE_DIR, top_n=G_S['top_n_p'] * G_S['top_n_s'], flag="",
            end_s_idx=G_S['end_s_idx'], species_path=G_S['species_path'])
    else:
        job_drivers.symbolic_path_2_real_path(
            FILE_DIR, top_n=G_S['top_n_p'], flag="",
            end_s_idx=None, species_path=G_S['species_path'])

    # copy SOHR/C++ routine files
    job_drivers.copy_sohr_files(FILE_DIR, species_path=G_S['species_path'])

    # # # species count
    # # job_drivers.species_count(FILE_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # reaction count
    # # job_drivers.reaction_count(FILE_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # initiation reaction count
    # # job_drivers.initiation_reaction_count(FILE_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species cycle
    # # job_drivers.species_cycle(FILE_DIR, top_n=G_S['top_n_p'], norm=True)

    # # # species production path
    # # job_drivers.species_production_path(
    # #     FILE_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # # # species production reaction
    # # job_drivers.species_production_reaction(
    # #     FILE_DIR, spe='OH', top_n=G_S['top_n_p'], norm=True)

    # propane make figures
    # job_drivers.propane_make_figures(
    #    FILE_DIR, species_path=G_S['species_path'])

    # send email
    job_drivers.send_email(FILE_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
