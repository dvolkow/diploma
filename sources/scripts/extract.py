#! /usr/bin/env python3

import sys
from astropy.io import fits


def get_idx_and_rs(filename):
    with fits.open(filename) as hdul:
        print(hdul.info())
        data = hdul[1].data
        idxs = data['APOGEE_ID'][:]
        rs   = data['RC_DIST'][:]
        return (idxs, rs)

def dump_cols(col1, col2, outpfname):
    f = open(outpfname, "w")
    for i in range(len(col1)):
        f.write(str(col1[i]) + ", " + str(col2[i]) + "\n")
    f.close()

if __name__ == "__main__":
    for param in sys.argv:
        print(param)
    idx, r = get_idx_and_rs(sys.argv[1])
    dump_cols(idx, r, sys.argv[2])



