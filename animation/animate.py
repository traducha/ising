#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
# import matplotlib
# matplotlib.use('Agg')
import os
import igraph as ig
import time
from numpy.random import random as rand
import numpy as np
import csv
from collections import OrderedDict
from pprint import pprint
from matplotlib import pyplot as plt
from main import *


name = 'test'
N = 50
M = 150
T = 20000


def plot_therm(n, m, time_step, mag, inc, largest_c, k_max):
    mag, = plt.plot(time_step, np.array(mag, dtype=float) / n, label='m')
    inc, = plt.plot(time_step, np.array(inc, dtype=float) / m, label='inc')
    largest_c, = plt.plot(time_step, np.array(largest_c, dtype=float) / n, label='S')
    k_max, = plt.plot(time_step, np.array(k_max, dtype=float) / n, label='k_max')
    handles = [mag, inc, largest_c, k_max]
    plt.legend(handles=handles, fontsize=7)
    plt.show()


def get_layout(g, niter=None, seed=None):
    return g.layout_graphopt(niter=niter, seed=seed, node_charge=0.0001, node_mass=50, spring_length=1)


def plot_g(g, plot_counter, layout, niter=None):
    ig.plot(g, 'plots/{}.png'.format(plot_counter), layout=layout, bbox=(720, 720))
    return get_layout(g, niter=niter, seed=layout), plot_counter + 1


def animate(g, n, m, times, layout, plot_counter=1):
    for j in xrange(50):
        layout, plot_counter = plot_g(g, plot_counter, layout, niter=5)

    mag, inc, largest_c, k_max, time_step = [], [], [], [], []
    count_changed = 0

    for i in xrange(times):
        g, changed = one_step(g, n, m)
        if changed:
            if count_changed % 5 == 0:
                print i
                mag.append(np.sum(g.vs()["spin"]))
                inc.append(count_incompatible(g))
                largest_c.append(len(max(g.clusters(), key=lambda clust: len(clust))))
                k_max.append(max(g.degree()))
                time_step.append(i)
                for j in xrange(5):
                    layout, plot_counter = plot_g(g, plot_counter, layout, niter=5)
            count_changed += 1

    plot_therm(n, m, time_step, mag, inc, largest_c, k_max)

    for j in xrange(200):
        layout, plot_counter = plot_g(g, plot_counter, layout, niter=5)

    return g


g = initialize_graph(N, M)
init_layout = get_layout(g, niter=500)
animate(g, N, M, T, init_layout)

os.system('ffmpeg -r 30 -pattern_type sequence -s 720x720 -start_number 0 -i "plots/%d.png" -q:v 1 {}.mp4'.format(name))
os.system('rm plots/*.png')
