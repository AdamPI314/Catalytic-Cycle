"""
main driver
"""

import os
import sys
import time

import update_settings as us
import job_drivers
import concatenate_results as cr
import read_write_configuration as rwc
import parse_spe_reaction_info as psri

if __name__ == '__main__':
    TIME_I = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    # TIME_Iter_counter = 0
    # TIME_Iter_N = 50
    # dt = 5.0e-9

    # while TIME_Iter_counter < TIME_Iter_N:
    #     # update settings
    #     us.update_settings(FILE_DIR, TIME_Iter_counter, dt=dt)

    #     # run jobs
    #     job_drivers.delete_srode_temp_files(FILE_DIR)
    #     try:
    #         job_drivers.path_ode_run(FILE_DIR)
    #     except RuntimeError:
    #         print("RuntimeError")
    #         exit()

    #     p_iter_N = rwc.read_iteration_Number(os.path.join(FILE_DIR, "input", "setting.json"))
    #     # organize files.  concatenate into one file
    #     cr.concatenate_time(FILE_DIR, p_iter_N=p_iter_N)
    #     cr.concatenate_concentration(FILE_DIR, p_iter_N=p_iter_N)
    #     cr.concatenate_temperature(FILE_DIR, p_iter_N=p_iter_N)
    #     cr.concatenate_pressure(FILE_DIR, p_iter_N=p_iter_N)

    #     TIME_Iter_counter += 1

    # run dlosde
    # job_drivers.run_dlsode(FILE_DIR)

    # make figures
    # job_drivers.make_figures(FILE_DIR)

    # send email
    # job_drivers.send_email(FILE_DIR)

    TIME_E = time.time()
    print("running time:\t" +
          str("{:.2f}".format((TIME_E - TIME_I) / 3600.0)) + " hours\n")
