#! /usr/bin/env python3
import math
import numpy as np
from statistics import stdev
import scipy.integrate as integrate


def my_erf(z):
    return math.sqrt(2. / math.pi) * integrate.quad(lambda x: math.pow(math.e, -0.5 * x ** 2), 0, z)[0]


'''
# 8
'''
def get_k(n):
    step = 0.00001
    current = step
    ideal = 1 - 1. / n
    while True:
        if (math.fabs(ideal - my_erf(current)) < step / 2):
            break;
        current += step
    return current



def distance_module(r):
    return 5 * math.log10(r)

def get_distance_scale(ds_1, ds_2):
    doubled = ds_1.viewkeys() & ds_2.viewkeys()
    N = len(doubled)
    delta_map = {}
    for idx in doubled:
        delta_map[idx] = distance_module(ds_1[idx]) - distance_module(ds_2[idx])
        
    mean_delta = reduce(lambda x, value:x + value, delta_map.values(), 0)
    mean_delta /= N
        
    mean_sq_delta = reduce(lambda x, value:x + value ** 2, delta_map.values(), 0)
    mean_sq_delta /= N

    sigma_delta = math.sqrt((mean_sq_delta - mean_delta ** 2) * N / (N - 1))

    sigma_delta_d = sigma_delta / math.sqrt(N)

    mean_r1r2 = math.pow(10, 0.2 * mean_delta)

    sigma_mean_r1r2 = math.log(10, math.e) / 5 * mean_r1r2 * sigma_delta_d 
    return (N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2, delta_map)

'''
# 9
'''
def get_ejected(k, delta_map, mean_delta, sigma_delta):
    res = {}
    for key, value in delta_map.items():
        if (math.fabs(value - mean_delta) > k * sigma_delta):
            res[key] = value

    return res
