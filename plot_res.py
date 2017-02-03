import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp


N = 100
L = 200
J = '1.000000'
h = '0.000000'
FI = '0.500000'

quants = ['mag', 'mag_abs', 'energy', 'incompatible', 'largest_clust', 'clust_num', 'largest_degree']

for q in quants:
	f_name = "res/{}_vs_B_N{}_L{}_J{}_h{}_FI{}.csv".format(q, N, L, J, h, FI)
	with open(f_name, 'rb') as file:
		x, value, std = csv.reader(file, delimiter=';', quotechar='|')

	for i in xrange(1, len(x)):
		x[i] = float(x[i])
		value[i] = float(value[i])
		std[i] = float(std[i])

	plt.errorbar(x[1:], value[1:], yerr=std[1:], fmt='o', fillstyle='none')
	plt.title("N={}, L={}, c={}".format(N, L, L/N))
	plt.xlabel(x[0])
	plt.ylabel(value[0])
	plt.savefig(f_name.replace(".csv" , ".png"), format="png")
	plt.clf()
