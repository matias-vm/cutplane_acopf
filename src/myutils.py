###############################################################################
#
# This code was written by Daniel Bienstock and is being maintained 
# by Matias Villagra.
# 
# Please report any bugs or issues (for sure there will be) to
#                         mjv2153@columbia.edu
#
# Oct 2023
###############################################################################


import sys


def breakexit(foo):
    stuff = input("("+foo+") break> ")
    if stuff == 'x' or stuff == 'q':
        sys.exit("bye")


