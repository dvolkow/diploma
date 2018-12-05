#! /usr/bin/env python3
import csv
import numpy as np
from math import sqrt 
import cds_core

EJECTED_F = "ejected.dat"
DELTA_F = "delta.dat"

def read_dataset(filename):
    f = open(filename, 'r')
    reader = csv.reader(f)
    res = {}
    for row in reader:
        res[str(row[0])] = float(row[1])

    print("Data has been readed, size = " + 
            str(len(res.keys())))
    return res


def dump_pairs(keys, filename, ds_1, ds_2):
    f = open(filename, "w")
    for key in keys:
        f.write("%f %f\n"% (ds_1[key], ds_2[key]))
    f.close()

def dump_res(N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2):
    f = open("cds_final_result.txt", "w")
    f.write("N               = %d\n"% (N))
    f.write("mean_delta      = %f\n"% (mean_delta))
    f.write("sigma_delta_d   = %f\n"% (sigma_delta_d))
    f.write("sigma_delta     = %f\n"% (sigma_delta))
    f.write("mean_r1r2       = %f\n"% (mean_r1r2))
    f.write("sigma_mean_r1r2 = %f\n"% (sigma_mean_r1r2))
    f.write("sgm             = %f\n"% (sigma_mean_r1r2 * sqrt(N)))
    f.close()
