import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint
from mean_field_2a import mean_field_phi as phi_mf
from mean_field_gamma import mean_field_gamma as gamma_mf
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'


n_list = [500]
G = '1.500000'
A = '0.000000'

quants = ['mag_abs', 'largest_degree', 'energy', 'incompatible', 'mag', 'largest_clust', 'clust_num', 'degree_corr']
names = [r'$|m|$', r'$k_{max}$', '$E$', 'incompatible', 'm', 'S', r'n_c', 'degree_corr']

os.chdir("../test_res")
g_mean = gamma_mf(500.0, 1500.0, 1.5, np.linspace(0.1, 15.0, 200))

plt.figure(figsize=(9, 8))
for j, q in enumerate(quants):
    if q == 'mag_abs':
        ax1 = plt.subplot(311)
        ax = ax1
    elif q == 'largest_degree':
        ax2 = plt.subplot(312, sharex=ax1)
        ax = ax2
    elif q == 'energy':
        ax3 = plt.subplot(313, sharex=ax2)
        ax = ax3
    else:
        continue

    for w, n in enumerate(n_list):
        m = n * 3
        norm = [1.0*n, 1.0*n, 1.0*m*n, 1.0*m, 1.0*n, 1.0*n, 1.0*n, 1.0]
        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, n, m, G, A)
        alpha, gamma = float(A), float(G)
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i])
            std[i] = float(std[i])

        ax.errorbar(x[1:], np.array(value[1:]) / norm[j], fmt='os^'[w], fillstyle='none')

    if q in ['largest_degree']:
        ax.set_ylim([-0.1, 1.1])
        ax.plot(np.linspace(0.1, 15.0, 200), g_mean[1], color='black', linewidth=2)
    if q in ['energy']:
        # ax.set_ylim([-1.0, 0.1])
        ax.set_xlabel(r'$T = 1/ \beta$', fontsize=16)
        ax.plot(np.linspace(0.1, 15.0, 200), g_mean[0], color='black', linewidth=2)
    if q in ['mag_abs']:
        ax.set_ylim([-0.01, 1.0])

    ax.set_ylabel(names[j], fontsize=16)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.tight_layout()
plt.subplots_adjust(hspace=0.000)
plt.show()
# plt.savefig('/home/tomasz/Desktop/1D_gamma.pdf', format="pdf")
plt.clf()
