import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint


N = 1000
M = 3000
G = '1.000000'
A = '0.000000'
temp_lim = [0, 150]
y_lim = [1, 3]
y_values = np.linspace(y_lim[0], y_lim[1], 81)[:-1]
y_var = 'G'  # 'A' or 'G'
aspect = temp_lim[1] / (y_lim[1] - y_lim[0])

quants = ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree', 'degree_corr']

os.chdir("../res/N1000/gamma")


for q in quants:
    # if q != 'mag_abs':
    #     continue
    value_matrix = []
    t_matrix = []
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

        # if g == '2.437500':
        #     break

        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, N, M, g, a)
        # print f_name
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i])
            std[i] = float(std[i])

        value_matrix.append(value[1:])
        t_matrix.append(x[1:])
        std_matrix.append(std[1:])

    # index = 8
    # final = [i[index] for i in value_matrix]
    # print t_matrix[5][index]
    # plt.plot(y_values, final)
    # plt.xlabel('$\phi$')
    # plt.ylabel('$|m|$')
    # plt.title('For T={}'.format(t_matrix[5][index]))
    # plt.show()
    # plt.savefig('/home/tomaszraducha/Pulpit/mag_vs_phi_T{}.png'.format(t_matrix[5][index]))

    im = plt.imshow(value_matrix, cmap=None, origin='lower', extent=temp_lim+y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.title('{} for $N={}, M={}$'.format(q, N, M))
    plt.xlabel('$T$')
    plt.ylabel('$\phi$' if y_var == 'A' else '$\gamma$')
    # plt.show()
    plt.savefig('/home/tomaszraducha/Pulpit/2D_{}_{}_{}{}.png'.format('PHI' if y_var == 'A' else 'GAMMA', q,
                                           'GA' if y_var == 'A' else 'AL', G if y_var == 'A' else A), format="png")
    plt.clf()

    im = plt.imshow(std_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none', aspect=aspect)
    plt.colorbar(im)
    plt.title('std of {} for $N={}, M={}$'.format(q, N, M))
    plt.xlabel('$T$')
    plt.ylabel('$\phi$' if y_var == 'A' else '$\gamma$')
    # plt.show()
    plt.savefig('/home/tomaszraducha/Pulpit/2D_STD_{}_{}_{}{}.png'.format('PHI' if y_var == 'A' else 'GAMMA', q,
                                               'GA' if y_var == 'A' else 'AL', G if y_var == 'A' else A), format="png")
    plt.clf()
