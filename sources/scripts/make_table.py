#! /usr/bin/env python3
import sys


class TableLine
    def __init__(self):
        self.R_0 = []
        self.A   = []
        self.u   = []
        self.v   = []
        self.w   = []
#       self.omega   = []
        self.theta   = []
        self.N   = 0


def read_dataset(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:

    return lines

def dump_res(filename, lines_array)
    with open(filename, 'w') as f:
#        f.write("$R_0$ & $" + lines_array[0][0] + "$ \\\hline")

data = read_dataset(sys.argv[1])
