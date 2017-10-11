import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp


N = 100
L = 300
J = '1.000000'
h = '0.000000'
G = '1.000000'
A = '1.000000'
quants = ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree', 'degree_corr']

os.chdir("../res")


def fetch_fi():
    fis = []
    pattern = re.compile(
        r'[a-zA-Z_]*_vs_B_N[0-9]{1,5}_L[0-9]{1,5}_J[0-9]\.[0-9]{6,6}_h0.000000_FI0.500000_GA([0-9]\.[0-9]{6,6})_AL([0-9]\.[0-9]{6,6})\.csv')
    for _file in glob.glob("*.csv"):
        print(_file)
        match = pattern.match(_file)
        fi = match.groups()[0]
        fis.append(fi)
    return list(set(fis))


for q in quants:
    f_name = "{}_vs_B_N{}_L{}_J{}_h{}_FI{}_GA{}_AL{}.csv".format(q, N, L, J, h, '0.500000', G, A)
    with open(f_name, 'rb') as file:
        x, value, std = csv.reader(file, delimiter=';', quotechar='|')

    for i in xrange(1, len(x)):
        x[i] = float(x[i])
        value[i] = float(value[i])
        std[i] = float(std[i])

    plt.errorbar(x[1:], value[1:], yerr=std[1:], fmt='o', fillstyle='none')
    plt.title("N={}, L={}, c={}".format(N, L, L / N))
    plt.xlabel(x[0])
    plt.ylabel(value[0])
    # plt.show()
    plt.savefig(f_name.replace(".csv", ".png"), format="png")
    plt.clf()
