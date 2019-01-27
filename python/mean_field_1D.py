# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
import scipy.special as sp
import scipy.stats as st

# params

N = 100.0
M = 300.0
gamma = 2.0
temp = np.linspace(0.1,60.0,1000)
beta = 1.0/temp
c = float(int(M/N))
energy = np.zeros(len(temp))
degree = np.zeros(len(temp))

exp_lim = 700.0

# functions

def L(i,N):
    return(i*(2.0*N-i-1.0)/2.0)

def k_a(i,N,M):
    return(i + 2.0*(M-L(i,N))/(N-i))

def log_A(i,N,M):
    An = (N - i)*(N - i - 1.0)/2.0
    Ak = M - L(i,N)
    return(sp.gammaln(An + 1.0) - sp.gammaln(Ak + 1.0) - sp.gammaln(An - Ak + 1.0))

def log_B(i,N,c):
    return(sp.gammaln(c + 1.0) - sp.gammaln(i + 1.0) + sp.gammaln(N - c + 1.0) - sp.gammaln(N - i + 1.0))

def log_exp(i,N,c,beta,gamma):
    return(beta*((i - c)*np.power(N-1.0,gamma) + (N-i)*np.power(k_a(i,N,M),gamma)))

def divide_exps(logs1,logs2):
    norm = np.max(logs2) # the highest log in the denominator
    counter = logs1 - norm
    denominator = logs2 - norm
    return(np.sum(np.exp(counter))/np.sum(np.exp(denominator)))

# calculations

km = np.repeat(N-1.0,c+1)
km[0] = st.poisson.ppf(1.0-1/N,2.0*M/N)
it = np.arange(c+1)
for b in range(len(beta)):
    logs_sum = np.zeros(c+1)
    logs_sum = log_B(it,N,c) + log_A(it,N,M) - log_A(c,N,M) + log_exp(it,N,c,beta[b],gamma)
    logs_sum_e = logs_sum + np.log(it) + gamma*np.log(N - 1.0) # energy counter
    logs_sum_e2 = logs_sum + np.log(N - it) + gamma*np.log(k_a(it,N,M)) # second energy counter
    logs_sum_d = logs_sum + np.log(km) # degree counter
    
    energy[b] = -(divide_exps(logs_sum_e,logs_sum) + divide_exps(logs_sum_e2,logs_sum))
    degree[b] = divide_exps(logs_sum_d,logs_sum)
energy = energy/(N*M)
degree = degree/(N-1.0)

# plots

py.figure(figsize=(12,9))
py.plot(temp,energy,'ro')
py.plot(temp,degree,'bo')



