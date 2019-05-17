#! /usr/bin/env python3
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stat
from math import sqrt
import matplotlib.pyplot as plt
import sys


FILE_OBJ=sys.argv[1]
obj = pd.read_csv(FILE_OBJ, delimiter = " ", names = ['R', 'Theta'])
plt.figure(figsize=(10, 10))

grid = sns.JointGrid(x = obj['R'], y = obj['Theta'], space = 0, height = 10, ratio = 8, xlim = [2, 15], ylim = [-400, 700])
grid.plot_joint(plt.scatter, color = "black", s = 0.1)
_ = grid.ax_marg_y.hist(obj['Theta'], color = "black", orientation = "horizontal", bins = np.arange(-400, 700, 10))
_ = grid.ax_marg_x.hist(obj['R'], color = "black", bins = np.arange(0, 15, 0.2))
FILE_LINE = sys.argv[2]
lines = pd.read_csv(FILE_LINE, delimiter = " ", names = ['R', 'Theta', 'LTheta', 'HTheta'])
sns.lineplot(x = lines['R'], y = lines['Theta'], palette = "tab10", linewidth = 2, color = 'b')
sns.lineplot(x = lines['R'], y = lines['LTheta'], palette = "tab10", linewidth = 0.5, color = 'b')
sns.lineplot(x = lines['R'], y = lines['HTheta'], palette = "tab10", linewidth = 0.5, color = 'b')

sublines = lines[['R', 'Theta', 'LTheta', 'HTheta']][lines['R'] > 6.5]
sublines = sublines[['R', 'Theta', 'LTheta', 'HTheta']][sublines['R'] < 10.6]
sns.lineplot(x = sublines['R'], y = sublines['Theta'], palette = "tab12", linewidth = 3, color = 'r')
max_p = sublines.loc[sublines['Theta'] == sublines['Theta'].max()]
min_p = sublines.loc[sublines['Theta'] == sublines['Theta'].min()]
print("min", min_p)
print("man", max_p)
sd = max_p['HTheta'].iloc[0] - max_p['LTheta'].iloc[0]
sd2 = min_p['HTheta'].iloc[0] - min_p['LTheta'].iloc[0]
delta = max_p['Theta'].iloc[0] - min_p['Theta'].iloc[0]


print("sqrt_sq = ", str(sqrt(sd ** 2 + sd2 ** 2)), ", fabs = ", abs(delta))

FILE_SUN = sys.argv[3]
sun = pd.read_csv(FILE_SUN, delimiter = " ", names = ['R', 'Theta'])
sns.scatterplot(x = sun['R'], y = sun['Theta'], color = 'yellow', size = 4)


plt.savefig("theta_obj.png")
