import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint
from mean_field import mean_field
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'


n_list = [500, 1000, 2000, 4000]
G = '1.000000'
A = '0.700000'

quants = ['energy', 'mag_abs', 'largest_degree', 'incompatible', 'mag', 'largest_clust', 'clust_num', 'degree_corr']
names = ['$E$', r'$|m|$', r'$\langle k \rangle$', 'incompatible', 'm', 'S', r'n_c', 'degree_corr']

os.chdir("../res_scaling/phi")

for j, q in enumerate(quants):
    plt.figure(figsize=(9, 3))
    plt.subplot(111)
    for n in n_list:
        m = n * 3
        norm = [1.0*m*n, 1.0*n, 1.0*n, 1.0*m, 1.0*n, 1.0*n, 1.0*n, 1.0]
        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, n, m, G, A)
        alpha, gamma = float(A), float(G)
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i])
            std[i] = float(std[i])

        plt.errorbar(x[1:], np.array(value[1:]) / norm[j], fmt='o')#, fillstyle='none')

    plt.xlabel(r'$T = 1/ \beta$', fontsize=16)
    plt.ylabel(names[j], fontsize=16)
    plt.tight_layout()
    # plt.show()
    plt.savefig('/home/tomasz/Desktop/scaling/phi/a_{}.pdf'.format(q), format="pdf")
    plt.clf()
