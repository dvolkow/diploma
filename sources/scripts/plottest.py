#! /usr/bin/env python3

import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stat
import matplotlib.pyplot as plt
import sys

# Below code for plotting ellipsods marginal distributions:
#pt = pd.read_table('../../build/result_1_10000_apogee-ucac_vr/R0Theta0.txt', delimiter=" ", names = ['x', 'y'], dtype = float)
#print(pt)
#g = (sns.jointplot(x = 'x', y = 'y', data = pt, kind = "kde", color="blue")).plot_joint(plt.scatter, c="b", s=1, linewidth=1, marker='+')
#g.set_axis_labels("$R_0$", "$\Theta$")
#plt.savefig("out.png")
#print(stat.pearsonr(pt['x'], pt['y']))


# View general and partial solutions:
PATH='../../build/test/'

def view_partial_sequences():
    i = str(sys.argv[2])
    vr_pt = pd.read_csv(PATH + i + "/vr_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    b_pt = pd.read_csv(PATH + i + "/b_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    l_pt = pd.read_csv(PATH + i + "/l_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    u_pt = pd.read_csv(PATH + i + "/u_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])

    sol = pd.concat([vr_pt, b_pt, l_pt, u_pt], axis = 0, sort=False, ignore_index=True)
    print(sol)


def view_united_sequences(size):
    sol = pd.DataFrame({'R_0': []})
    for j in range(1,size):
        u_pt = pd.read_csv(PATH + str(j) + "/u_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
        sol = pd.concat([sol, u_pt], axis = 0, sort=False, ignore_index=True)
    print(sol)

def view_residuals():
    pt = pd.read_csv('../../build/test/residuals.txt', delimiter = " ", names = ['V_R', 'mu_b', 'mu_l'])
    g = sns.distplot(pt['V_R'], kde = True, rug = True, bins = 100)
    plt.savefig("VR.png")
    plt.clf()
    g = sns.distplot(pt['mu_b'], kde = True, rug = True, bins = 100)
    plt.savefig("mu_b.png")
    plt.clf()
    g = sns.distplot(pt['mu_l'], kde = True, rug = True, bins = 100)
    plt.savefig("mu_l.png")

if sys.argv[1] == "s":
    view_partial_sequences()
elif sys.argv[1] == "v":
    view_united_sequences(int(sys.argv[2]))
elif sys.argv[1] == "r":
    view_residuals()
else: 
    print("Nothing. Keys: s, v, r.")





