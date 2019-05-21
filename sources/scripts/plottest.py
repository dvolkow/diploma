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
    #print(pt)
    g = (sns.jointplot(x = 'x', y = 'y', data = pt, kind = "kde", color="blue")).plot_joint(plt.scatter, c="b", s=1, linewidth=1, marker='+')
    g.set_axis_labels("$R_0$", "$\Theta$")
    plt.savefig("out.png")
    print(stat.pearsonr(pt['x'], pt['y']))


def view_errors():
    PATH=str(sys.argv[2])
    sol = pd.DataFrame({})
    for j in range(1, int(sys.argv[3])):
        pt = pd.read_csv(PATH + str(j) + "/" + str(sys.argv[4]), delimiter = " ")
        sol = pd.concat([sol, pt], axis = 0, sort = False, ignore_index = True)
    print(sol)



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
        u_pt = pd.read_csv(PATH + str(j) + "/" + str(sys.argv[4]), delimiter=" ", names = ['N', 'R_0', 'SD', 'u', 'v', 'w', 'omega_0', 'A', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6', 'theta7', 'theta8', 'theta9', 'theta10'])
        sol = pd.concat([sol, u_pt], axis = 0, sort=False, ignore_index=True)
    print(sol)

def view_cormatrix():
    FILE=str(sys.argv[2])
    pt = pd.read_csv(FILE, delimiter = " ")
    #print(pt)
    g = sns.PairGrid(pt)
    g.map_diag(sns.kdeplot)
    g.map_lower(sns.kdeplot, n_levels = 6)
    g.map_upper(plt.scatter, s = 0.5)
    plt.savefig("pairplot.png")

def view_unf():
    FILE=str(sys.argv[2])
    na = pd.read_csv(FILE, delimiter = " ")
    print(na.T.to_latex())

def view_r0res():
    FILE=str(sys.argv[2])
    na = pd.read_csv(FILE, delimiter = " ")
    for index, row in na.iterrows():
        print("$_{-" + str(row[0]) + "}^{+" + str(row[1]) + "}$")



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

    pt = pd.read_csv(PATH + '/theta_errs.txt', delimiter = " ", names = ['R', 'dT'])
    g = sns.distplot(pt['dT'], kde = True, rug = True, bins = 100)
    plt.savefig("theta_err.png")
    plt.clf()

    pt = pd.read_csv(PATH + '/uni_profile.txt', delimiter = " ", names = ['R_0', 'Sigma'])
    #print(pt)
    g = sns.lineplot(x = "R_0", y="Sigma", data=pt, markers=True, linewidth=1.5, palette="tab10")
    plt.savefig("profile.png")

def view_profile():
    plt.figure(figsize=(10, 10))
#    sns.JointGrid(x = obj['R'], y = obj['Theta'], space = 0, height = 10, ratio = 8, xlim = [2, 15], ylim = [-400, 700])
    for i in range(1,int(sys.argv[3])):
        pt = pd.read_csv(str(sys.argv[2]) + str(i) + str(sys.argv[4]), delimiter = " ", names = ['R_0', 'Sigma'])
        sns.lineplot(x = "R_0", y="Sigma", data=pt, markers=True, linewidth=0.6, palette="tab10")
    plt.savefig("profiles.png")

def view_lb():
    fig = plt.figure()
    plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot()
    pt = pd.read_csv(str(sys.argv[2]), delimiter = " ", names = ['ra', 'dec', 'l', 'b', 'nv', 'vr', 'vr_err', 'r', 'pm_ra', 'pm_dec', 'pm_ra_err', 'pm_dec_err', 'match'])
    ax1 = (sns.jointplot(x = 'l', y = 'b', data = pt, xlim = [0, 360], ratio = 5, ylim = [-90, 90], s=0.1, color="black")).plot_joint(plt.scatter, c="gray", s=0.1, linewidth=1, marker='+')
    ax1.set_axis_labels("l, deg", "b, deg")
    plt.savefig("lb.png")
    plt.clf()
    ax1 = sns.distplot(pt['vr'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("vr_distr.png")
    plt.clf()
    ax1 = sns.distplot(pt['vr_err'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("vr_err_distr.png")
    plt.clf()
    ax1 = sns.distplot(pt['pm_ra'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("pm_ra_distr.png")
    plt.clf()
    ax1 = sns.distplot(pt['pm_ra_err'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("pm_ra_err_distr.png")
    plt.clf()
    ax1 = sns.distplot(pt['pm_dec'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("pm_dec_distr.png")
    plt.clf()
    ax1 = sns.distplot(pt['pm_dec_err'], kde = True, rug = False, bins = 100, color = 'black')
    plt.savefig("pm_dec_err_distr.png")
    plt.clf()


def view_xyz():
    PATH=sys.argv[2]
#    NAME = sys.argv[3]
#    NAME_ERR = sys.argv[4]
    fig = plt.figure()
    plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot()

    pt = pd.read_csv(PATH + '/final_xyz.txt', delimiter = " ", names = ['X', 'Y', 'Z'])

    ax1 = (sns.jointplot(x = 'X', y = 'Y', data = pt, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="black")).plot_joint(plt.scatter, c="gray", s=0.1, linewidth=1, marker='+')
    ax1.set_axis_labels("X, kpc", "Y, kpc")
    plt.savefig("XYobj.png")

    ax1 = (sns.jointplot(x = 'X', y = 'Z', data = pt, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="black")).plot_joint(plt.scatter, c="gray", s=0.1, linewidth=1, marker='+')
    ax1.set_axis_labels("X, kpc", "Z, kpc")
    plt.savefig("XZobj.png")

#    ax1 = (sns.jointplot(x = 'Y', y = 'Z', data = pt, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="b")).plot_joint(plt.scatter, c="b", s=0.1, linewidth=1, marker='+')
#    ax1.set_axis_labels("Y", "Z")
#    plt.savefig("YZobj.png")

#    perrors = pd.read_csv(NAME, delimiter = " ", names = ['X', 'Y', 'Z'])

#    ax2 = fig.add_subplot()

#    ax2 = (sns.jointplot(x = 'X', y = 'Y', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="green")).plot_joint(plt.scatter, c="green", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("X", "Y")
#    plt.savefig("XYmiss.png")

#    ax2 = (sns.jointplot(x = 'X', y = 'Z', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="green")).plot_joint(plt.scatter, c="green", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("X", "Z")
#    plt.savefig("XZmiss.png")

#    ax2 = (sns.jointplot(x = 'Y', y = 'Z', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="green")).plot_joint(plt.scatter, c="green", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("Y", "Z")
#    plt.savefig("YZmiss.png")

#    perrors = pd.read_csv(NAME_ERR, delimiter = " ", names = ['X', 'Y', 'Z'])
#    ax2 = (sns.jointplot(x = 'X', y = 'Y', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="red")).plot_joint(plt.scatter, c="red", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("X", "Y")
#    plt.savefig("XYerr.png")

#    ax2 = (sns.jointplot(x = 'X', y = 'Z', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="red")).plot_joint(plt.scatter, c="red", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("X", "Z")
#    plt.savefig("XZerr.png")

#    ax2 = (sns.jointplot(x = 'Y', y = 'Z', data = perrors, xlim = [-10, 10], ylim = [-10, 10], s=0.1, color="red")).plot_joint(plt.scatter, c="red", s=0.1, linewidth=1, marker='+')
#    ax2.set_axis_labels("Y", "Z")
#    plt.savefig("YZerr.png")



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
elif sys.argv[1] == "th":
    view_one_pairplot()
elif sys.argv[1] == "err":
    view_errors()
elif sys.argv[1] == "p":
    view_profile()
elif sys.argv[1] == "a":
    view_unf()
elif sys.argv[1] == "r0":
    view_r0res()
elif sys.argv[1] == "lb":
    view_lb()
else: 
    print("Nothing. Keys: s, v, r.")





