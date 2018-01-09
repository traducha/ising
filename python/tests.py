# -*- coding: utf-8 -*-
from matplotlib import rc
rc('font', **{'family' :'serif'})
import matplotlib.pyplot as plt
from scipy.special import gamma as g
from scipy.special import binom
from scipy.special import hyp2f1
import numpy as np
import math
import os
import csv
from mean_field import mean_field
from mean_field_phi import mean_field_phi


N = 100.0
M = 400.0
t_max = 60

gammas = np.linspace(0.0, 3.0, 100)
temp = np.linspace(0.1, t_max, 100)

trans_temp = []
first = True

for g in []:#gammas:
    _, degrees = mean_field(N, M, g, temp)
    t = 0.0
    again = True
    once = False
    for i, degree in enumerate(degrees):
        if degree > 55 and again:
            once = True
            t = temp[i]
            if i == len(degrees) - 1:
                if first is False:
                    t = 1000
                else:
                    first = False
        elif once:
            again = False
    trans_temp.append(t)

final_g, final_t = [], []
for i, t in enumerate(trans_temp):
    if 0.0 < t < 1000:
        final_t.append(t)
        final_g.append(gammas[i])

# plt.plot(final_t, final_g, color='#FFFF33', linewidth=2)

phis = np.linspace(0.0, 2.0, 1000)
temp = np.linspace(0.1, t_max, 1000)

trans_temp = []
first = True

for p in phis:
    _, degrees = mean_field_phi(N, M, p, temp)
    t = 0.0
    for i, degree in enumerate(degrees):
        if degree < 18:
            t = temp[i]
            break
    # print p, t
    # print degrees
    trans_temp.append(t)

final_p, final_t = [], []
for i, t in enumerate(trans_temp):
    if 0.0 < t < 1000:
        final_t.append(t)
        final_p.append(phis[i])

plt.plot(final_t, final_p, color='#FFFF33', linewidth=2)

plt.plot([0, t_max], [1.477, 1.477], color='#006600', linewidth=4)

G = '1.000000'
A = '0.000000'
temp_lim = [0, t_max]
y_lim = [0, 2]
y_values = np.arange(101) * float(y_lim[1]) / 101.0
y_var = 'A'  # 'A' or 'G'
aspect = temp_lim[1] / y_lim[1]

os.chdir("../res_c{}/res_phi".format(int(M/N)))
quant = 'largest_degree'  # ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree', 'degree_corr']

if 1:
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

        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(quant, int(N), int(M), g, a)
        # print f_name
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
    # plt.title('{} for $N={}, M={}$'.format(quant, N, M))
    plt.xlabel('$T$', fontsize=16)
    plt.ylabel('$\phi$' if y_var == 'A' else '$\gamma$', fontsize=16)
    # plt.show()
    plt.savefig('/home/tomasz/Desktop/2D_{}_{}_{}{}.pdf'.format('PHI' if y_var == 'A' else 'GAMMA', quant,
                                           'GA' if y_var == 'A' else 'AL', G if y_var == 'A' else A), format="pdf")
    plt.clf()




