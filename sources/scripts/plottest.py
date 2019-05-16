#! /usr/bin/env python3

import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stat
import matplotlib.pyplot as plt
import sys

# Below code for plotting ellipsods marginal distributions:
def view_one_pairplot():
    FILE=str(sys.argv[2])
    pt = pd.read_table(FILE, delimiter=" ", names = ['x', 'y'], dtype = float)
    print(pt)
    g = (sns.jointplot(x = 'x', y = 'y', data = pt, kind = "kde", color="blue")).plot_joint(plt.scatter, c="b", s=1, linewidth=1, marker='+')
    g.set_axis_labels("$R_0$", "$\Theta$")
    plt.savefig("out.png")
    print(stat.pearsonr(pt['x'], pt['y']))




# View general and partial solutions:

def view_partial_sequences():
    PATH=str(sys.argv[3])
    i = str(sys.argv[2])
    vr_pt = pd.read_csv(PATH + i + "/vr_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    b_pt = pd.read_csv(PATH + i + "/b_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    l_pt = pd.read_csv(PATH + i + "/l_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
    u_pt = pd.read_csv(PATH + i + "/u_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])

    sol = pd.concat([vr_pt, b_pt, l_pt, u_pt], axis = 0, sort=False, ignore_index=True)
    print(sol)


def view_united_sequences(size):
    PATH=str(sys.argv[3])
    sol = pd.DataFrame({'R_0': []})
    for j in range(1,size):
        u_pt = pd.read_csv(PATH + str(j) + "/u_unfresult.txt", delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
        sol = pd.concat([sol, u_pt], axis = 0, sort=False, ignore_index=True)
    print(sol)

def view_cormatrix():
    FILE=str(sys.argv[2])
    pt = pd.read_csv(FILE, delimiter = " ")
    print(pt)
    g = sns.PairGrid(pt)
    g.map_diag(sns.kdeplot)
    g.map_lower(sns.kdeplot, n_levels = 6)
    g.map_upper(plt.scatter, s = 0.5)
    plt.savefig("pairplot.png")

def view_residuals():
    PATH=sys.argv[2]
    pt = pd.read_csv(PATH + '/residuals.txt', delimiter = " ", names = ['V_R', 'mu_b', 'mu_l'])
    g = sns.distplot(pt['V_R'], kde = True, rug = True, bins = 100)
    plt.savefig("VR.png")
    plt.clf()
    g = sns.distplot(pt['mu_b'], kde = True, rug = True, bins = 100)
    plt.savefig("mu_b.png")
    plt.clf()
    g = sns.distplot(pt['mu_l'], kde = True, rug = True, bins = 100)
    plt.savefig("mu_l.png")

    plt.clf()
    pt = pd.read_csv(PATH + '/uni_profile.txt', delimiter = " ", names = ['R_0', 'Sigma'])
    print(pt)
    g = sns.lineplot(x = "R_0", y="Sigma", data=pt, markers=True, linewidth=1.5, palette="tab10")
    plt.savefig("profile.png")

def view_xyz():
    PATH=sys.argv[2]
    fig = plt.figure()
    ax1 = fig.add_subplot()
    pt = pd.read_table(PATH + '/final_xyz.txt', delimiter = " ", names = ['X', 'Y', 'Z'])
    ax1 = (sns.jointplot(x = 'X', y = 'Y', data = pt, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="b")).plot_joint(plt.scatter, c="b", s=0.1, linewidth=1, marker='+')
    ax1.set_axis_labels("X", "Y")
    perrors = pd.read_table(PATH + '/ERROR_LIMITED', delimiter = " ", names = ['X', 'Y', 'Z'])
    plt.savefig("XYobj.png")
    ax2 = fig.add_subplot()
    ax2 = (sns.jointplot(x = 'X', y = 'Y', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="red")).plot_joint(plt.scatter, c="red", s=0.1, linewidth=1, marker='+')
    ax2.set_axis_labels("X", "Y")
    plt.savefig("XY.png")

if sys.argv[1] == "s":
    view_partial_sequences()
elif sys.argv[1] == "v":
    view_united_sequences(int(sys.argv[2]))
elif sys.argv[1] == "r":
    view_residuals()
elif sys.argv[1] == "mk":
    view_cormatrix()
elif sys.argv[1] == "xy":
    view_xyz()
else: 
    print("Nothing. Keys: s, v, r.")





