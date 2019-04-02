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


name = 'ising'
params = {
    'n': 50,
    'm': 150,
    'temp': None,
    'beta': None,
    'gamma': 0.0,
    'phi': 0.0,
    'times': None,
    'agg': None,
}
params_list = [
    dict(params, agg=10, times=500, temp=0.5, beta=1.0/0.5),
    dict(params, agg=10, times=8000, temp=0.5, beta=1.0/0.5, gamma=2.0),
    dict(params, agg=10, times=1000, temp=50.0, beta=1.0/50.0, gamma=2.0),
    dict(params, agg=45, times=3500, temp=50.0, beta=1.0/50.0, gamma=3.0),
    dict(params, agg=30, times=5000, temp=5.0, beta=1.0/5.0, phi=1.0),
    dict(params, agg=40, times=6000, temp=5.0, beta=1.0/5.0, phi=0.5),
    dict(params, agg=20, times=1200, temp=50.0, beta=1.0/50.0, phi=0.5),
    dict(params, agg=15, times=12000, temp=50.0, beta=1.0/50.0, phi=2.0),
]


def plot_therm(n, m, time_step, mag, inc, largest_c, k_max, params=None):
    mag, = plt.plot(time_step, np.array(mag, dtype=float) / n, label='m')
    inc, = plt.plot(time_step, np.array(inc, dtype=float) / m, label='inc')
    largest_c, = plt.plot(time_step, np.array(largest_c, dtype=float) / n, label='S')
    k_max, = plt.plot(time_step, np.array(k_max, dtype=float) / n, label='k_max')
    handles = [mag, inc, largest_c, k_max]
    plt.legend(handles=handles, fontsize=7)
    plt.show()


def get_layout(g, niter=None, seed=None, params=None):
    return g.layout_graphopt(niter=niter, seed=seed, node_charge=0.0001, node_mass=50, spring_length=1)


def plot_g(g, plot_counter, layout, niter=None, params=None):
    ig.plot(g, 'plots/{}.png'.format(plot_counter), layout=layout, bbox=(720, 720))
    return get_layout(g, niter=niter, seed=layout), plot_counter + 1


def animate(g, n, m, times, layout, plot_counter=1, agg=None, params=None):
    for j in xrange(50):
        layout, plot_counter = plot_g(g, plot_counter, layout, niter=5)

    mag, inc, largest_c, k_max, time_step = [], [], [], [], []
    count_changed = 0

    for i in xrange(times):
        g, changed = one_step(g, n, m, params=params)
        if changed:
            if count_changed % agg == 0:
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

    for j in xrange(150):
        layout, plot_counter = plot_g(g, plot_counter, layout, niter=5)

    return g, layout, plot_counter


def video(g, layout, plot_counter, name, divided=False, params=None):
    start_number = plot_counter
    g, layout, plot_counter = animate(g, params['n'], params['m'], params['times'], layout, plot_counter=plot_counter,
                                      agg=params['agg'], params=params)
    if divided:
        os.system('ffmpeg -r 30 -pattern_type sequence -s 720x720 -start_number {} '
                  '-i "plots/%d.png" -c:v libx264 -q:v 1 vids/{}.mp4'.format(start_number, name))
    return g, layout, plot_counter


if __name__ == '__main__':
    g = initialize_graph(params['n'], params['m'], params=params)
    layout = get_layout(g, niter=500, params=params)
    plot_counter = 1

    for i, params in enumerate(params_list):
        g, layout, plot_counter = video(g, layout, plot_counter, ''.join([name, str(i+1)]), divided=True, params=params)

    # os.system('ffmpeg -r 30 -pattern_type sequence -s 720x720 -start_number 0 -i "plots/%d.png" -c:v libx264 -q:v 1 vids/{}3.mp4'.format(name))
    os.system('rm plots/*.png')
