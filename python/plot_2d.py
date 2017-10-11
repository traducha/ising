import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint


N = 100
M = 300
G = '2.000000'
A = '1.000000'
y_values = np.arange(100) / 50.0
y_var = 'A'  # 'A' or 'G'
temp_lim = [0, 120]
y_lim = [0, 2]
aspect = temp_lim[1] / y_lim[1]

quants = ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree', 'degree_corr']

os.chdir("../res_phi2")


for q in quants:
    value_matrix = []
    std_matrix = []

    for y in y_values:
        g = G
        a = A
        if y_var == 'A':
            a = str(y).ljust(8, '0')
        elif y_var == 'G':
            g = str(y).ljust(8, '0')
        else:
            raise

        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, N, M, g, a)
        print f_name
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i])
            std[i] = float(std[i])

        value_matrix.append(value[1:])
        std_matrix.append(std[1:])

    im = plt.imshow(value_matrix, cmap=None, origin='lower', extent=temp_lim+y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.title('{} for $N={}, M={}$'.format(q, N, M))
    plt.xlabel('$T$')
    plt.ylabel('$\phi$' if y_var == 'A' else '$\gamma$')
    # plt.show()
    plt.savefig('2D_{}_{}_{}{}.png'.format('PHI' if y_var == 'A' else 'GAMMA', q,
                                           'GA' if y_var == 'A' else 'AL', G if y_var == 'A' else A), format="png")
    plt.clf()

    im = plt.imshow(std_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.title('std of {} for $N={}, M={}$'.format(q, N, M))
    plt.xlabel('$T$')
    plt.ylabel('$\phi$' if y_var == 'A' else '$\gamma$')
    # plt.show()
    plt.savefig('2D_STD_{}_{}_{}{}.png'.format('PHI' if y_var == 'A' else 'GAMMA', q,
                                               'GA' if y_var == 'A' else 'AL', G if y_var == 'A' else A), format="png")
    plt.clf()
