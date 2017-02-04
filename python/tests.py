# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.special import gamma as g
from scipy.special import binom
from scipy.special import hyp2f1
import numpy as np
import math


N = 20.0
L = 40.0
J = 1.0
h = 0.0


def Zn(B, n):
    m = (n - (N - n)) / N
    exp_ = np.exp(L * J * B + N * h * m * B)
    
    au = n * (N - n)
    ae = 0.5 * (n * (n - 1) + (N - n) * (N - n - 1))
    frac = (math.factorial(ae)) / (math.factorial(L) * math.factorial(ae - L))

    fi = hyp2f1(-au, -L, 1 + ae - L, np.exp(-2 * J * B))
    
    return exp_ * frac * fi


def Z(B):
    return np.sum([binom(N, i) * Zn(B, i) for i in xrange(0, int(N+1))])


def Zn_po_h(B, n):
    m = (n - (N - n)) / N
    return Zn(B, n) *  N * m * B


def Z_po_h(B):
    return np.sum([binom(N, i) * Zn_po_h(B, i) for i in xrange(0, int(N+1))])


def mag(B):
    return Z_po_h(B) / (B * Z(B))


temp = np.linspace(0.1, 10.0, 10)
magnet = [mag(1.0/t) * 10000000000000000 for t in temp]  # 10000000000000000 żeby wykres lepiej wyglądał

print temp
print magnet

plt.scatter(temp, magnet)
plt.xlabel('Temp.')
plt.ylabel('Mag.')
plt.show()
plt.clf()














