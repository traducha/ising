# -*- coding: utf-8 -*-
import csv
import numpy as np
from math import factorial as fac
import scipy.special as sp
import matplotlib.pyplot as plt

# parametry

#B = 1. # mam nadzieje, ze to nie absurdalny wybor

N = 20.0
L = N*(N-1)/2.0
J = 1.
h = 0.
c = L/N

def E(n):
    m = (2.*n - N)/N
    return -J * (L - 2 * n * (N - n)) - N * h * m

def Z1(B):
    return np.sum([(fac(N) / (fac(n) * fac(N-n))) * np.exp(-B*E(n)) for n in xrange(0,int(N)+1)])

def Z_po_h1(B):
    return np.sum([(fac(N) / (fac(n) * fac(N-n))) * np.exp(-B*E(n)) * B*N*(np.abs((2.*n - N)/N)) for n in xrange(0,int(N)+1)])

def mag1(B):
    return Z_po_h1(B) / (B * N * Z1(B))


temp = np.linspace(0.1, 40.0, 200)
magnet = [mag1(1.0/t) for t in temp]  # 10000000000000000 żeby wykres lepiej wyglądał

plt.figure(figsize=(12,8))
plt.scatter(temp, magnet)
#plt.scatter(temp, magnet2,color='red')
#plt.scatter(temp, magnet3,color='green')
plt.xlabel('Temp.')
plt.ylabel('Mag.')

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
N = 20.0
L = 90.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

# obliczenia
N = 20.0
L = 80.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

# obliczenia
N = 20.0
L = 60.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

# obliczenia
N = 20.0
L = 40.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

# obliczenia
N = 20.0
L = 20.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

# obliczenia
N = 20.0
L = 10.0

temp = np.linspace(0.1, 40.0, 200)
magnet = [np.abs(mag(1.0/t)) for t in temp]
plt.scatter(temp, magnet, color='r')

if 1:
		f_name = "res/mag_vs_B_N10L140.csv"
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i]) / N
		    std[i] = float(std[i]) / N

		plt.errorbar(x[1:], value[1:], yerr=std[1:], fmt='o', fillstyle='none')
		plt.title("N={}, L={}, c={}".format(N, L, L/N))

plt.show()
plt.clf()
