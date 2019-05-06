#! /usr/bin/env python3
import math
import numpy as np
from statistics import stdev
from functools import reduce
import scipy.integrate as integrate
from scipy.optimize import fsolve

G_EJECTED_OBJS = []

def my_erf(z):
    return math.sqrt(2. / math.pi) * integrate.quad(lambda x: math.pow(math.e, -0.5 * x ** 2), 0, z)[0]

def funct(x, n, p):
    return p - (1 - my_erf(x)) * n

'''
# 8, 10B
@p = 0.05 or 1
'''
def get_k(n, p):
    res = fsolve(func = funct, x0 = [0], args = (n, p), xtol = 1e-8)
    return res[0]


def distance_module(r):
    return 5 * math.log10(r)

def get_distance_scale(ds_1, ds_2):
    doubled = ds_1.keys() & ds_2.keys()
    delta_map = {}
    for idx in doubled:
        delta_map[idx] = distance_module(ds_1[idx]) - distance_module(ds_2[idx])
    return delta_map

def solution(delta_map):
    N = len(delta_map)
    mean_delta = reduce(lambda x, value:x + value, delta_map.values(), 0)
    mean_delta /= N
        
    mean_sq_delta = reduce(lambda x, value:x + value ** 2, delta_map.values(), 0)
    mean_sq_delta /= N
    # print("mean_sq_delta = ", mean_sq_delta)

    sigma_delta = math.sqrt((mean_sq_delta - mean_delta ** 2) * N / (N - 1))
    # print("sigma_delta = ", sigma_delta)

    sigma_delta_d = sigma_delta / math.sqrt(N)
    # print("sigma_delta_d = ", sigma_delta_d)

    mean_r1r2 = math.pow(10, 0.2 * mean_delta)
    # print("mean_r1r2 = ", mean_r1r2)

    sigma_mean_r1r2 = math.log(10, math.e) / 5 * mean_r1r2 * sigma_delta_d 
    # print("sigma_mean_r1r2 = ", sigma_mean_r1r2)
    return (N, mean_delta, sigma_delta_d, sigma_delta, mean_r1r2, sigma_mean_r1r2, delta_map)

'''
# 9
'''
def get_ejected(k, delta_map, mean_delta, sigma_delta):
    res = {}
    for key, value in delta_map.items():
        if (math.fabs(value - mean_delta) > k * sigma_delta):
            res[key] = math.fabs(value)
    return res

def max_item_key(ejected):
    if len(ejected) == 0:
        return None
    finded = False
    res_v = 0
    res_k = 0
    for key, value in ejected.items():
        if not finded or res_v < value:
            res_k, res_v = key, value
            finded = True
    return res_k


def get_new_delta_map(L, ejected, delta_map, sigma_delta):
    ej_len = len(ejected)
    if ej_len == 0:
        return delta_map
    tmp_ej = {}

    N = len(delta_map)
    for i in range(L - 1):
        key_max = max_item_key(ejected)
        # print("finded key_max = ", key_max)
        tmp_ej[key_max] = ejected[key_max]
        ejected.pop(key_max)
        ej_len -= 1
        if ej_len == 0:
            break
    for key in ejected.keys():
        delta_map.pop(key)
        G_EJECTED_OBJS.append(key)

    key_max = max_item_key(tmp_ej)
    if key_max is None:
        return delta_map
    if (tmp_ej[key_max] > get_k(N, 0.05) * sigma_delta):
        delta_map.pop(key_max)
        G_EJECTED_OBJS.append(key)

    return delta_map


