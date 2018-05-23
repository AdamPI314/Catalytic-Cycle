"""
rename files with regular expression
"""

import os
import sys
import time
import fnmatch


def parse_filename(data_dir, re_expr, fmt="png", prefix="C3H8_H_"):
    """
    parse filename with regular expression
    """
    for _, _, files in os.walk(data_dir):
        for f_n in files:
            if fnmatch.fnmatch(f_n, "*" + re_expr + "*." + fmt):
                os.rename(os.path.join(data_dir, f_n),
                          os.path.join(data_dir, prefix + f_n))

    return


if __name__ == '__main__':
    INIT_TIME = time.time()

    DATA_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir, os.pardir, "SOHR_DATA"))

    # RE_EXPR = "straight_"
    # FMT = "png"
    # RE_EXPR = "network_H"
    # FMT = "gexf"
    RE_EXPR = "network_H"
    FMT = "gephi"
    PREFIX = "C3H8_H_"
    parse_filename(os.path.join(DATA_DIR, "output"), RE_EXPR, FMT, PREFIX)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
