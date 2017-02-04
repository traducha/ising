# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

# parametry

#B = 1. # mam nadzieje, ze to nie absurdalny wybor

N = 100.#500.
L = 200.#1000.
J = 1.
h = 0.09
c = L/N

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
    log_Z_sum = np.sum([np.exp(log_Z_n(B,n_u) - log_Z_0 + log_c_n(n_u) - log_c_0) for n_u in range(0,int(N+1))])
    log_Z_sum += np.exp(log_Z_n(B,N) - log_Z_0)
#    return(log_Z_0 + log_c_0 + np.log(1. + log_Z_sum))
    return(np.log(1. + log_Z_sum)) # wywalilem log_Z_0 i log_c_0, poniewaz sie skracaja

def Z_po_h(B):
    log_Z_0 = log_Z_n(B,0)
    log_c_0 = log_c_n(0)
    Z_po_h_sum = np.sum([(2.*n_u - N)*np.exp(log_Z_n(B,n_u) - log_Z_0 + log_c_n(n_u) - log_c_0)/N for n_u in range(0,int(N+1))])
    Z_po_h_sum += np.exp(log_Z_n(B,N) - log_Z_0) + 1.0
    return(N * B * Z_po_h_sum) # bez Z_0 i c_0, ktore by sie skrocily

def mag(B):
    return(Z_po_h(B) / (B * N * np.exp(log_Z(B))))

# obliczenia

temp = np.linspace(0.1, 10.0, 200)
magnet = [mag(1.0/t) for t in temp]  # 10000000000000000 żeby wykres lepiej wyglądał

#print(temp)
#print(magnet)

# obrazek

plt.figure(figsize=(12,8))
plt.scatter(temp, magnet)
#plt.scatter(temp, magnet2,color='red')
#plt.scatter(temp, magnet3,color='green')
plt.xlabel('Temp.')
plt.ylabel('Mag.')

#L = 300
#temp = np.linspace(0.1, 20.0, 200)
#magnet = [mag(1.0/t) for t in temp]  # 10000000000000000 żeby wykres lepiej wyglądał
#plt.scatter(temp, magnet)

#L = 600
#temp = np.linspace(0.1, 20.0, 200)
#magnet = [mag(1.0/t) for t in temp]  # 10000000000000000 żeby wykres lepiej wyglądał
#plt.scatter(temp, magnet)

plt.show()
plt.clf()

# kolejny krok to policzenie sumy...


# testowanie poprawnosci danych wejsciowych

#for n_u in range(int(N+1)):
#    a_e = (n_u*(n_u - 1) + (N - n_u)*(N - n_u - 1))/2.
#    a_u = n_u*(N - n_u)
#    if(a_e<L):
#        print('a_e: ' + str(n_u))
#    if(sp.hyp2f1(-a_u,-L,1.+a_e-L,np.exp(-2.*J*B))==np.inf):
#        print('2f1: ' + str(n_u))

