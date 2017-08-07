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

    END_TIME = 3.51e-4

    # # run dlosde
    # job_drivers.run_dlsode(FILE_DIR, END_TIME)

    # # run monte carlo trajectory
    # job_drivers.run_mc_trajectory(FILE_DIR, END_TIME)

    # # convert symbolic pathway to real pathway
    # # with real species names and real reaction expression
    # job_drivers.symbolic_path_2_real_path(FILE_DIR, top_n=50, flag="")

    # evaluate path integral-->pathway probability
    job_drivers.evaluate_pathway_probability(FILE_DIR, top_n=5, num=1, flag="")

    # # send email
    # job_drivers.send_email(FILE_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
