# -*- coding: utf-8 -*-
import numpy as np
import pylab as py
import scipy.special as sp
import scipy.stats as st

# params

N = 100.0
M = 300.0
phi = 1.0
temp = np.linspace(0.1, 40.0, 1000)
beta = 1.0 / temp
n_m = np.ceil(0.5 + np.sqrt(1.+8.*M)/2.)
av_k = 2.0 * M / N
energy = np.zeros(len(temp))
degree = np.zeros(len(temp))


# functions

def log_2(i, N):
    return (N - i + 1) * np.log(2)


def log_A(i, N):
    return sp.gammaln(N + 1.0) - sp.gammaln(i + 1.0) - sp.gammaln(N - i + 1.0)


def log_B(i, M):
    As = i * (i - 1.0) / 2.0
    return sp.gammaln(As + 1.0) - sp.gammaln(M + 1.0) - sp.gammaln(As - M + 1.0)


def log_exp(i, M, beta, phi):
    return beta * np.power(M, 2.0 * phi + 1.0) * np.power(2.0 / i, 2.0 * phi) / av_k


def log_E(i, M, phi):
    return (2.0 * phi + 1.0) * np.log(M) + 2.0 * phi * np.log(2.0 / i) - np.log(av_k)


# 1D calculations

def mean_field_phi(n, m, fi, temp):
    global N
    global M
    global beta
    global av_k
    global n_m
    global phi
    global energy
    global degree
    N = n
    M = m
    phi = fi
    beta = 1.0 / temp
    av_k = 2.0 * M / N
    energy = np.zeros(len(temp))
    degree = np.zeros(len(temp))
    n_m = np.ceil(0.5 + np.sqrt(1. + 8. * M) / 2.)

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
        energy[b] = -np.sum(np.exp(logs_sum_E)) / np.sum(np.exp(logs_sum))
        if (np.isnan(energy[b])):
            energy[b] = -np.exp(np.max(logs_sum_E) - np.max(logs_sum))
        # degree
        degree[b] = np.sum(np.exp(logs_sum_D)) / np.sum(np.exp(logs_sum))
        if (np.isnan(degree[b])):
            degree[b] = np.exp(np.max(logs_sum_D) - np.max(logs_sum))

    return energy, degree


if __name__ == '__main__':
    phi = 1.0
    temp = np.linspace(0.1, 400, 1000)
    energy, degree = mean_field_phi(1000.0, 3000.0, phi, temp)
    idx = np.argwhere(np.diff(np.sign(np.array([1 if x == 1 else -1 for x in np.array(degree < 18)]))) != 0).reshape(-1) + 0
    print temp[idx]

    # py.plot(temp, energy, 'r')
    # py.show()

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
        degree[p, b] = np.sum(np.exp(logs_sum_D)) / np.sum(np.exp(logs_sum))
        if (np.isnan(degree[p, b])):
            degree[p, b] = np.exp(np.max(logs_sum_D) - np.max(logs_sum))

# degree phase diagram plot

py.figure(figsize=(12, 12))
py.imshow(degree, cmap=None, origin='lower', extent=[0.1, 500.0, 0.0, 1.5], aspect=300)
py.colorbar()
py.xlabel('temperature')
py.ylabel('phi')
py.show()
"""