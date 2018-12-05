#! /usr/bin/env python3

import sys
import cds_core as core
import cds_io as io
#from cds_wrappers import cds_solution_show

if __name__ == '__main__':
    counter = 1
    refdata_1 = io.read_dataset(sys.argv[1])
    refdata_2 = io.read_dataset(sys.argv[2])
    N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2, delta_map = core.solution(core.get_distance_scale(refdata_1, refdata_2))
    old_size = len(delta_map)
    print("#%d: %f %f (size %d)"% (counter, mean_r1r2, sigma_mean_r1r2, old_size))
    ejected = core.get_ejected(core.get_k(N, 1), delta_map, mean_delta, sigma_delta)
    delta_map = core.get_new_delta_map(1, ejected, delta_map, sigma_delta)
    new_size = len(delta_map)
    counter += 1
    while new_size != old_size:
        N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2, delta_map = core.solution(delta_map)
        print("#%d: %f %f (size %d)"% (counter, mean_r1r2, sigma_mean_r1r2, new_size))
        counter += 1
        old_size = new_size
        ejected = core.get_ejected(core.get_k(N, 1), delta_map, mean_delta, sigma_delta)
        delta_map = core.get_new_delta_map(1, ejected, delta_map, sigma_delta)
        new_size = len(delta_map)
    print("\nFINAL RESULTS: ")
    print("mean_delta      = %f"% (mean_delta))
    print("sigma_delta_d   = %f"% (sigma_delta_d))
    print("sigma_delta     = %f"% (sigma_delta))
    print("mean_r1r2       = %f"% (mean_r1r2))
    print("sigma_mean_r1r2 = %f"% (sigma_mean_r1r2))
    io.dump_pairs(core.G_EJECTED_OBJS, io.EJECTED_F, refdata_1, refdata_2)
    io.dump_pairs(delta_map.keys(), io.DELTA_F, refdata_1, refdata_2)
    io.dump_res(N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2)
