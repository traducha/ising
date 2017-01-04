import matplotlib.pyplot as plt
import sys
import csv
import math
from scipy.special import gamma as g
import numpy as np

N = 5.0
L = 10.0
J = 1.0
h = 0.0



temp = np.linspace(1.0, 10.0, 100)

t = temp[0]
B = 1.0/t


def Zn(B, n):
    m = 1.0 * (n - (N - n)) / N
    exp_ = np.exp(L * J * B + N * h * m * B)
    
    ae = 0.5 * (n * (n-1) + (N-n) * (N-n-1))
    print ae
    frac = (g(1+ae)) / (g(1+L) * g(1 + ae -L))

    return 1.0


def Z(B):
    res = 0.0
    for i in xrange(0, int(N+1)):
        res += Zn(B, i)
    return res


print Z(B)
