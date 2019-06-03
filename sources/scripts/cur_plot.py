#! /usr/bin/env python3
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stat
from math import sqrt
import matplotlib.pyplot as plt
import sys


FILE_OBJ=sys.argv[1]
fig = plt.figure()
obj = pd.read_csv(FILE_OBJ, delimiter = " ", names = ['R', 'Theta'])
plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot()

grid = sns.JointGrid(x = obj['R'], y = obj['Theta'], space = 0, height = 10, ratio = 8, xlim = [0, 18], ylim = [-150, 500])
grid.plot_joint(plt.scatter, color = "black", s = 0.1)
_ = grid.ax_marg_y.hist(obj['Theta'], color = "black", orientation = "horizontal", bins = np.arange(-150, 500, 10))
_ = grid.ax_marg_x.hist(obj['R'], color = "black", bins = np.arange(0, 18, 0.2))
FILE_LINE = sys.argv[2]
lines = pd.read_csv(FILE_LINE, delimiter = " ", names = ['R', 'Theta', 'LTheta', 'HTheta'])
sns.lineplot(x = lines['R'], y = lines['Theta'], palette = "tab10", linewidth = 2, color = 'b')
sns.lineplot(x = lines['R'], y = lines['LTheta'], palette = "tab10", linewidth = 0.5, color = 'b')
sns.lineplot(x = lines['R'], y = lines['HTheta'], palette = "tab10", linewidth = 0.5, color = 'b')

low_R = 5
high_R = 12
low_th = 0
high_th = 400

sublines = lines[['R', 'Theta', 'LTheta', 'HTheta']][lines['R'] > low_R]
sublines = sublines[['R', 'Theta', 'LTheta', 'HTheta']][sublines['R'] < high_R]
#sns.lineplot(x = sublines['R'], y = sublines['Theta'], palette = "tab12", linewidth = 3, color = 'r')
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
#ax1.set_axis_label(['R, kpc', r'$\Theta$'])
plt.savefig('theta_obj.png')

plt.clf()
grid = sns.JointGrid(x = obj['R'], y = obj['Theta'], space = 0, height = 10, ratio = 15, xlim = [low_R, high_R], ylim = [low_th, high_th])
grid.plot_joint(plt.scatter, color = "black", s = 0.1)
#_ = grid.ax_marg_y.hist(obj['Theta'], color = "black", orientation = "horizontal", bins = np.arange(low_th, high_th, 5))
_ = grid.ax_marg_x.hist(obj['R'], color = "black", bins = np.arange(low_R, high_R, 0.1))
FILE_LINE = sys.argv[2]
lines = pd.read_csv(FILE_LINE, delimiter = " ", names = ['R', 'Theta', 'LTheta', 'HTheta'])
sns.lineplot(x = lines['R'], y = lines['Theta'], palette = "tab10", linewidth = 2, color = 'b')
sns.lineplot(x = lines['R'], y = lines['LTheta'], palette = "tab10", linewidth = 0.5, color = 'b')
sns.lineplot(x = lines['R'], y = lines['HTheta'], palette = "tab10", linewidth = 0.5, color = 'b')
sns.scatterplot(x = sun['R'], y = sun['Theta'], color = 'yellow', size = 10)
plt.savefig("zoomed.png")


