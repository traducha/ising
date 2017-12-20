import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'

N = 1000
M = 3000
G = '1.000000'
A = '0.000000'
temp_lim = [0, 150]
y_lim = [1, 3]
y_values = np.linspace(y_lim[0], y_lim[1], 81)[:-1]
y_var = 'G'  # 'A' or 'G'
aspect = temp_lim[1] / (y_lim[1] - y_lim[0])

quants = ['energy', 'mag_abs', 'largest_degree', 'incompatible', 'mag', 'largest_clust', 'clust_num', 'degree_corr']
norm = [1.0*N*M, 1.0*N, 1.0*N, 1.0*M, 1.0*N, 1.0*N, 1.0*N, 1.0]

os.chdir("../res_2D/gamma")


for j, q in enumerate(quants):
    value_matrix = []
    std_matrix = []

    for y in y_values:
        g = G
        a = A
        if y_var == 'A':
            a = str(round(y, 6)).ljust(8, '0')
        elif y_var == 'G':
            g = str(round(y, 6)).ljust(8, '0')
        else:
            raise

        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, N, M, g, a)
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i]) / norm[j]
            std[i] = float(std[i]) / norm[j]

        value_matrix.append(value[1:])
        std_matrix.append(std[1:])

    im = plt.imshow(value_matrix, cmap=None, origin='lower', extent=temp_lim+y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.xlabel(r'$T = 1/ \beta$', fontsize=16)
    plt.ylabel(r'$\phi$' if y_var == 'A' else r'$\gamma$', fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig('/home/tomasz/Desktop/2D/gamma/g_2D{}.pdf'.format(q), format='pdf')
    plt.clf()

    im = plt.imshow(std_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.xlabel(r'$T = 1/ \beta$', fontsize=16)
    plt.ylabel(r'$\phi$' if y_var == 'A' else r'$\gamma$', fontsize=20)
    plt.tight_layout()
    # plt.show()
    plt.savefig('/home/tomasz/Desktop/2D/gamma/g_STD_2D{}.pdf'.format(q), format='pdf')
    plt.clf()
