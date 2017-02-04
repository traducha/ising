import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp


N = 100.0
L = 200.0
J = 1.0
h = 0.0

# Funkcje

def log_Z_n(B,n_u):
    m = (2.*n_u - N)/N
    a_e = (n_u*(n_u - 1.) + (N - n_u)*(N - n_u - 1.)) * 0.5
    a_u = n_u*(N - n_u)
    
    exp_ln = (L*J*B + N*h*m*B)
    fact_ln = sp.gammaln(1.+a_e)-sp.gammaln(1.+L)-sp.gammaln(1.+a_e-L)
    hyper_ln = np.log(sp.hyp2f1(-a_u,-L,1.+a_e-L,np.exp(-2.*J*B)))
    
    return(exp_ln + fact_ln + hyper_ln)

def log_c_n(n_u):
    return(sp.gammaln(N + 1) - sp.gammaln(N - n_u + 1) - sp.gammaln(n_u + 1))

def log_Z(B):
    log_Z_0 = log_Z_n(B,0)
    log_c_0 = log_c_n(0)
    log_Z_sum = np.sum([np.exp(log_Z_n(B,n_u) - log_Z_0 + log_c_n(n_u) - log_c_0) for n_u in range(1,int(N))])
    log_Z_sum += np.exp(log_Z_n(B,N) - log_Z_0)
#    return(log_Z_0 + log_c_0 + np.log(1. + log_Z_sum))
    return(np.log(1. + log_Z_sum)) # wywalilem log_Z_0 i log_c_0, poniewaz sie skracaja

def Z_po_h(B):
    log_Z_0 = log_Z_n(B,0)
    log_c_0 = log_c_n(0)
    Z_po_h_sum = np.sum([np.abs(2.*n_u - N)*np.exp(log_Z_n(B,n_u) - log_Z_0 + log_c_n(n_u) - log_c_0)/N for n_u in range(1,int(N))])
    Z_po_h_sum += np.exp(log_Z_n(B,N) - log_Z_0) + 1.0
    return(N * B * Z_po_h_sum) # bez Z_0 i c_0, ktore by sie skrocily

def mag(B):
    return(Z_po_h(B) / (B * N * np.exp(log_Z(B))))

# obliczenia

temp = np.linspace(0.1, 10.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]


if 1:
		f_name = "res/mag_vs_B_N100L200.csv"
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i]) / N
		    std[i] = float(std[i]) / N

		plt.errorbar(x[1:], value[1:], yerr=std[1:], fmt='o', fillstyle='none')
		plt.title("N={}, L={}, c={}".format(N, L, L/N))
		plt.xlabel(x[0])
		plt.ylabel(value[0])


plt.plot(temp, magnet, color='black')

plt.show()
#plt.savefig("algorithms_comparsion_mag.png", format="png")
plt.clf()
