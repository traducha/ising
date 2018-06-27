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


n_list = [500, 750, 1000]
N = 1000.0
N2 = 500.0
N3 = 750.0
G = '1.500000'
A = '0.000000'
mean_field_lim = (0.1, 12.0, 200)

quants = ['mag_abs', 'largest_degree', 'energy', 'incompatible', 'mag', 'largest_clust', 'clust_num', 'degree_corr']
names = [r'$|m|$', r'$k_{max}$', '$E$', 'incompatible', 'm', 'S', r'n_c', 'degree_corr']

os.chdir("../res_scaling_fin/gamma")
g_mean = gamma_mf(float(N), N*3.0, float(G), np.linspace(*mean_field_lim))
g_mean_2 = gamma_mf(float(N2), N2*3.0, float(G), np.linspace(*mean_field_lim))
g_mean_3 = gamma_mf(float(N3), N3*3.0, float(G), np.linspace(*mean_field_lim))

plt.figure(figsize=(9 / 1.2, 8 / 1.2))
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
        norm = [1.0*n, 1.0*n, 1.0*n*m, 1.0*m, 1.0*n, 1.0*n, 1.0*n, 1.0]
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
        ax.plot(np.linspace(*mean_field_lim), g_mean[1], '--', color='red', linewidth=1)
        ax.plot(np.linspace(*mean_field_lim), g_mean_2[1], color='blue', linewidth=1)
        ax.plot(np.linspace(*mean_field_lim), g_mean_3[1], '-.', color='green', linewidth=1)
    if q in ['energy']:
        # ax.set_ylim([-1.0, 0.1])
        ax.set_xlabel(r'$T = 1/ \beta$', fontsize=18)
        ax.plot(np.linspace(*mean_field_lim), g_mean[0], '--', color='red', linewidth=1)
        ax.plot(np.linspace(*mean_field_lim), g_mean_2[0], color='blue', linewidth=1)
        ax.plot(np.linspace(*mean_field_lim), g_mean_3[0], '-.', color='green', linewidth=1)
    if q in ['mag_abs']:
        ax.set_ylim([-0.01, 1.0])

    ax.set_ylabel(names[j], fontsize=18)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.tight_layout()
plt.subplots_adjust(hspace=0.000)
# plt.show()
plt.savefig('/home/tomasz/Desktop/1D_gamma.pdf', format="pdf")
plt.clf()
