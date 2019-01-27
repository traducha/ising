# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
import scipy.special as sp
import scipy.stats as st

# params

N = 0
M = 0
phi = 0
temp = 0
beta = 0
n_m = 0
k_av = 0

energy = 0
degree = 0

exp_lim = 700.0


# functions

def log_2(i, N):
    return ((N - i + 1) * np.log(2))


def log_A(i, N):
    return (sp.gammaln(N + 1.0) - sp.gammaln(i + 1.0) - sp.gammaln(N - i + 1.0))


def log_B(i, M):
    As = i * (i - 1.0) / 2.0
    return (sp.gammaln(As + 1.0) - sp.gammaln(M + 1.0) - sp.gammaln(As - M + 1.0))


def log_exp(i, M, beta, phi):
    return beta * np.power(M, 2.0 * phi + 1.0) * np.power(2.0 / i, 2.0 * phi) * np.power(1.0 / k_av, phi)


def log_E(i, M, phi):
    return (2.0 * phi + 1.0) * np.log(M) + 2.0 * phi * np.log(2.0 / i) + phi * np.log(1.0 / k_av)


def divide_exps(logs1, logs2):
    norm = np.max(logs2)  # the highest log in the denominator
    counter = logs1 - norm
    denominator = logs2 - norm
    return (np.sum(np.exp(counter)) / np.sum(np.exp(denominator)))


# 1D calculations

def mean_field(n, m, ph, temps):
    global N
    global M
    global beta
    global temp
    global n_m
    global k_av
    global phi
    global energy
    global degree
    N = n
    M = m
    phi = ph
    temp = temps
    beta = 1.0 / temp
    k_av = 2.0 * M / N
    n_m = np.ceil(0.5 + np.sqrt(1. + 8. * M) / 2.)
    energy = np.zeros(len(temp))
    degree = np.zeros(len(temp))

    it = np.arange(n_m, N + 1)
    k_s = np.minimum(st.poisson.ppf(1.0 - 1 / it, 2.0 * M / it), it - 1)
    for b in range(len(beta)):
        logs_norm = log_2(n_m, N) + log_A(n_m, N) + log_B(n_m, M) + log_exp(n_m, M, beta[b],
                                                                            phi)  # log of the first element of Z sum
        # we subtract logs_norm which is equivalent of dividing both the counter and denominator by a constant factor
        # this way elements in both counter and denominator are much smaller but the result is the same (numerical trick)
        logs_sum = log_2(it, N) + log_A(it, N) + log_B(it, M) + log_exp(it, M, beta[b], phi) - logs_norm
        logs_sum_E = logs_sum + log_E(it, M, phi)  # counter for energy
        logs_sum_D = logs_sum + np.log(k_s)  # counter for degree
        # energy
        energy[b] = -divide_exps(logs_sum_E, logs_sum) / (M * N)
        # degree
        degree[b] = divide_exps(logs_sum_D, logs_sum) / N
    return energy, degree

# 1D plots

if __name__ == '__main__':
    temperature = np.linspace(0.1, 70.0, 1000)
    energy, degree = mean_field(500.0, 1500.0, 0.7, temperature)
    py.figure(figsize=(9, 6))
    py.plot(temp, energy, 'ro')
    py.show()

    py.figure(figsize=(9, 6))
    py.plot(temp, degree, 'bo')
    py.show()

# degree phase diagram
"""
temp = np.linspace(0.1, 500.0, 1000)
beta = 1.0 / temp
phis = np.linspace(0.0, 1.5, 300)

degree = np.zeros((len(phis), len(temp)))

it = np.arange(n_m, N + 1)
k_s = np.minimum(st.poisson.ppf(1.0 - 1 / it, 2.0 * M / it), it - 1)
for b in range(len(beta)):
    for p in range(len(phis)):
        logs_norm = log_2(n_m, N) + log_A(n_m, N) + log_B(n_m, M) + log_exp(n_m, M, beta[b], phis[
            p])  # log of the first element of Z sum
        # we subtract logs_norm which is equivalent of dividing both the counter and denominator by a constant factor
        # this way elements in both counter and denominator are much smaller but the result is the same (numerical trick)
        logs_sum = log_2(it, N) + log_A(it, N) + log_B(it, M) + log_exp(it, M, beta[b], phis[p]) - logs_norm
        logs_sum_D = logs_sum + np.log(k_s)  # counter for degree
        # degree
        degree[p, b] = divide_exps(logs_sum_D, logs_sum)

# degree phase diagram plot

py.figure(figsize=(12, 12))
py.imshow(degree, cmap=None, origin='lower', extent=[0.1, 500.0, 0.0, 1.5], aspect=300)
py.colorbar()
py.xlabel('temperature')
py.ylabel('phi')
py.show()
"""
