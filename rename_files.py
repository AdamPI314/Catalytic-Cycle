"""
rename files with regular expression
"""

import os
import sys
import time


def parse_filename(file_dir, re_expr, format="png"):
    """
    parse filename with regular expression
    """
    print(file_dir)

    return

if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(
        sys.argv[0]), os.pardir, os.pardir, os.pardir))
    print(FILE_DIR)

    RE_EXPR = "straight_H"
    FORMAT = "png"
    parse_filename(os.path.join(FILE_DIR, "output"), RE_EXPR, FORMAT)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
