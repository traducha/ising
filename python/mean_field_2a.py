# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
import scipy.special as sp
import scipy.stats as st

# params

N = 100.0
M = 300.0
phi = 0.7
temp_max = 40
temp = np.linspace(0.1, temp_max, 1000)
beta = 1.0 / temp
n_m = np.ceil(0.5 + np.sqrt(1. + 8. * M) / 2.)

energy = np.zeros(len(temp))
degree = np.zeros(len(temp))

exp_lim = 700.0


# functions

def log_2(i, N):
    return ((N - i + 1) * np.log(2))


def log_A(i, N):
    return (sp.gammaln(N + 1.0) - sp.gammaln(i + 1.0) - sp.gammaln(N - i + 1.0))


def log_B(i, M):
    As = i * (i - 1.0) / 2.0
    return (sp.gammaln(As + 1.0) - sp.gammaln(M + 1.0) - sp.gammaln(As - M + 1.0))


def log_exp(i, M, N, beta, phi):
    return (beta * np.power(M, 2.0 * phi + 1.0) * np.power(2.0 / i, 2.0 * phi) / np.power(2 * M / N, phi))


def log_E(i, M, N, phi):
    return ((2.0 * phi + 1.0) * np.log(M) + 2.0 * phi * np.log(2.0 / i) - phi * np.log(2 * M / N))


def divide_exps(logs1, logs2):
    norm = np.max(logs2)  # the highest log in the denominator
    counter = logs1 - norm
    denominator = logs2 - norm
    return (np.sum(np.exp(counter)) / np.sum(np.exp(denominator)))


# 1D calculations

def mean_field_phi(n, m, fi, temp):
    global N
    global M
    global beta
    global n_m
    global phi
    global energy
    global degree
    N = n
    M = m
    phi = fi
    beta = 1.0 / temp
    energy = np.zeros(len(temp))
    degree = np.zeros(len(temp))
    n_m = np.ceil(0.5 + np.sqrt(1. + 8. * M) / 2.)

    it = np.arange(n_m, N + 1)
    k_s = np.minimum(st.poisson.ppf(1.0 - 1.0 / it, 2.0 * M / it), it - 1)
    for b in range(len(beta)):
        logs_sum = log_2(it, N) + log_A(it, N) + log_B(it, M) + log_exp(it, M, N, beta[b], phi)
        logs_sum_E = logs_sum + log_E(it, M, N, phi)  # counter for energy
        logs_sum_D = logs_sum + np.log(k_s)  # counter for degree
        # energy
        energy[b] = -divide_exps(logs_sum_E, logs_sum) / (N * M)
        # degree
        degree[b] = divide_exps(logs_sum_D, logs_sum) / N

    return energy, degree


# 1D plots

if __name__ == '__main__':
    # py.figure(figsize=(9, 6))
    py.plot(temp, energy, 'ro')
    py.show()

    # py.figure(figsize=(9, 6))
    py.plot(temp, degree, 'bo')
    py.show()

