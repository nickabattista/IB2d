#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import sys
# Path Reference to where Driving code is found #
sys.path.append('IBM_Blackbox')
import IBM_Driver as Driver

def main2d(sim_path):
    """Function reads all input parameters and submits all information
    to the IBM_Driver files.
    """
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Frontend for pyIB2d. '
                                    'Use file-select dialog or terminal arg '
                                    'to point to simulation directory.')
    parser.add_argument('-d','--fdiag', action='store_true', help='Use file-dialog')
    parser.add_argument('-p', '--path', help='Path to input2d for sim')
    args = parser.parse_args()
    if args.path:
        print('Path option selected.')
        #main(args.path)
    elif args.fdiag:
        print('file-dialog option selected.')
        pass
    else:
        print('Invalid flags. See -h.')
