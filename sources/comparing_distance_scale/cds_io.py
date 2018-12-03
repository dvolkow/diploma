#! /usr/bin/env python3
import csv
import numpy as np
from math import radians
import cds_core

def read_dataset(filename):
    data = cds_core.Data()
    f = open(filename, 'r')
    reader = csv.reader(f)
    res = {}
    for row in reader:
        res[str(row[0])] = float(row[1])

    print("Data has been readed, size = " + 
            str(len(res.keys())))
    return res


