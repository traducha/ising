# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.special import gamma as g
from scipy.special import binom
from scipy.special import hyp2f1
import numpy as np
import math
import os
import csv


N = 100.0
M = 300.0
c = M / N
nh = np.floor(c)
n = 0.5 + np.sqrt(1.0 + 8.0 * M) / 2.0

def f1(x):
    return x * np.log((n - 1.0)**2.0 / (N - 1)) + np.log(N)


def f2(x):
    return np.log(((nh - 1.0) / 2.0) * (N - 1)**x + (N - nh) * nh**x)


approx = np.log((nh - 1.0) / (2.0 * N)) / (2.0 * np.log(n / N))
print 'approx:', approx

phis = np.linspace(0.0, 3.0, 1000)
f1_ = np.array([f1(p) for p in phis])
f2_ = np.array([f2(p) for p in phis])

idx = np.argwhere(np.diff(np.sign(f1_ - f2_)) != 0).reshape(-1) + 0
print 'exact', phis[idx]

plt.plot(phis, f1_, 'blue')
plt.plot(phis, f2_, 'green')
plt.plot(phis[idx], f1_[idx], 'ro')
plt.xlabel('phi')
plt.show()
