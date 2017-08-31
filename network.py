"""
network related, generate network, save to json, save to gephi compatible file
"""

import os
import time
import sys
import networkx


def init_network(file_dir):
    """
    init network
    """
    return

if __name__ == '__main__':
    INIT_TIME = time.time()

    FILE_DIR = os.path.abspath(os.path.join(os.path.realpath(sys.argv[0]), os.pardir, os.pardir))
    print(FILE_DIR)

    END_TIME = time.time()

    print("Time Elapsed:\t{:.5} seconds".format(END_TIME - INIT_TIME))
