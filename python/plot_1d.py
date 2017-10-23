import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint
from mean_field import mean_field


N = 100
M = 400
G = 1.0
A = 1.2

quants = ['energy']  # ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree', 'degree_corr']
norm = 1.0 * M * N

os.chdir("../res_c4/res_phi")
diff = 1000.0
alpha, gamma = 0.0, 0.0

for q in quants:
    pattern = re.compile(
        r'{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA([0-9]\.[0-9]*)_AL([0-9]\.[0-9]*)\.csv'.format(q, N, M))

    for _file in glob.glob("{}*.csv".format(q)):
        match = pattern.match(_file)
        ga, al = match.groups()

        new_diff = abs(A - float(al)) + abs(G - float(ga))
        if new_diff < diff:
            f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, N, M, ga, al)
            diff = new_diff
            alpha, gamma = float(al), float(ga)

    print f_name
    with open(f_name, 'rb') as file:
        x, value, std = csv.reader(file, delimiter=';', quotechar='|')

    for i in xrange(1, len(x)):
        x[i] = float(x[i])
        value[i] = float(value[i])
        std[i] = float(std[i])

    temp = np.linspace(x[1], x[-1], 1000)
    # mean_e, mean_d = mean_field(N, M, gamma, temp)

    plt.errorbar(x[1:], np.array(value[1:]) / norm, yerr=np.array(std[1:]) / norm, fmt='o', fillstyle='none')
    # plt.plot(temp, mean_e, 'black')
    plt.title("{}, N={}, M={}, ga={}, phi={}".format(q, N, M, gamma, alpha))
    plt.xlabel('temp.')
    plt.ylabel(q)
    plt.show()
    # plt.savefig(f_name.replace(".csv", ".png"), format="png")
    plt.clf()
