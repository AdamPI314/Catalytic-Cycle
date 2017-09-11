"""
main driver
"""

import os
import sys
import time

import job_drivers

if __name__ == '__main__':
    TIME_I = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    END_TIME = 0.800000099855441071
    TOP_N = 2000
    TAU = 0.5
    ATOM_FOLLOWED = "C"
    INIT_SPE = 62
    END_SPE = ""

    # # run dlosde
    # job_drivers.run_dlsode(FILE_DIR, END_TIME)

    # # write specie concentration at a time to file
    # job_drivers.spe_concentration_at_time_w2f(FILE_DIR, tau=TAU)

    # # run monte carlo trajectory
    # job_drivers.run_mc_trajectory(
    #     FILE_DIR, END_TIME, n_traj=10000000, atom_followed=ATOM_FOLLOWED, init_spe=INIT_SPE, max_tau=TAU)

    # # evaluate path integral-->pathway probability
    job_drivers.evaluate_pathway_probability(
        FILE_DIR, top_n=TOP_N, num=1, flag="", n_traj=10000, atom_followed=ATOM_FOLLOWED, init_spe=INIT_SPE, max_tau=TAU, top_s_n=10)

    # copy SOHR/C++ routine files
    job_drivers.copy_sohr_files(FILE_DIR)

    # # convert symbolic pathway to real pathway
    # # with real species names and real reaction expression
    # job_drivers.symbolic_path_2_real_path(
    #     FILE_DIR, top_n=TOP_N, flag="", end_spe=END_SPE)

    # # species count
    # job_drivers.species_count(FILE_DIR, top_n=TOP_N, norm=True)

    # # reaction count
    # job_drivers.reaction_count(FILE_DIR, top_n=TOP_N, norm=True)

    # # initiation reaction count
    # job_drivers.initiation_reaction_count(FILE_DIR, top_n=TOP_N, norm=True)

    # # species cycle
    # job_drivers.species_cycle(FILE_DIR, top_n=TOP_N, norm=True)

    # # species production path
    # job_drivers.species_production_path(
    #     FILE_DIR, spe='OH', top_n=TOP_N, norm=True)

    # # species production reaction
    # job_drivers.species_production_reaction(
    #     FILE_DIR, spe='OH', top_n=TOP_N, norm=True)

    # send email
    job_drivers.send_email(FILE_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
