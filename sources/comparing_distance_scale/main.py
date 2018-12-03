#! /usr/bin/env python3

import cds_core as core
import cds_io as io
from cds_wrappers import cds_solution_show

if __name__ == '__main__':
    refdata_1 = io.read_dataset()
    refdata_2 = io.read_dataset()
    res = core.get_distance_scale(refdata_1, refdata_2)
    cds_solution_show(res)
