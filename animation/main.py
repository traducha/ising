#!/usr/bin/python
# -*- coding: utf-8 -*-
import igraph as ig
import time
from numpy.random import random as rand
import numpy as np
import csv
from collections import OrderedDict
from matplotlib import pyplot as plt


GREEN = '#26A57C'
RED = '#DA3A49'
edge_width = 2

params = {
    'n': None,  # 100
    'm': None,  # 300
    'temp': None,  # 0.5
    'beta': None,  # 1.0/temp
    'gamma': None,  # 2.0
    'phi': None,  # 0.0
    'thermal': None,  # 300000
    'sample_time': None,  # 1000
}


def initialize_graph(n, m, params=None):
    g = ig.Graph.Erdos_Renyi(n=n, m=m)
    g.vs()["spin"] = (np.int_(rand((n, 1)) * 2) * 2) - 1
    g.vs()["color"] = [GREEN if spin == 1 else RED for spin in g.vs["spin"]]
    g.es()['width'] = edge_width
    return g


def spin(g, index, params=None):
    return g.vs(index)["spin"][0][0]


def energy_change_spin(g, spin_change, v_index, neigs, params=None):
    delta = 0.0
    v_degree = g.degree(v_index)
    neigs_degree = g.degree(neigs)

    for i in xrange(len(neigs)):
        delta += (neigs_degree[i] ** params['phi']) * g.vs(neigs[i])["spin"][0][0]
    delta *= - spin_change * (v_degree ** params['phi'])
    return delta


def energy_change_edge_alpha(g, edge_list, old_from, old_to, new_from, new_to, params=None):
    old_from_d, old_to_d, new_from_d, new_to_d = g.degree([old_from, old_to, new_from, new_to])
    delta = 0.0
    tmp_delta = 0.0

    # i,j -> m,n
    delta += spin(g, old_from) * spin(g, old_to) * (old_from_d ** params['phi']) * (old_to_d ** params['phi'])\
        - spin(g, new_from) * spin(g, new_to) * ((new_from_d + 1) ** params['phi']) * ((new_to_d + 1) ** params['phi'])

    # old_from neighbours
    neigs = g.neighbors(old_from)
    neigs_d = g.degree(neigs)
    for i, neig in enumerate(neigs):
        if neig not in (old_from, old_to, new_from, new_to):
            tmp_delta += (neigs_d[i] ** params['phi']) * spin(g, neig)
    delta += tmp_delta * spin(g, old_from) * (old_from_d ** params['phi'] - (old_from_d - 1) ** params['phi'])
    tmp_delta = 0.0

    # old_to neighbours
    neigs = g.neighbors(old_to)
    neigs_d = g.degree(neigs)
    for i, neig in enumerate(neigs):
        if neig not in (old_from, old_to, new_from, new_to):
            tmp_delta += (neigs_d[i] ** params['phi']) * spin(g, neig)
    delta += tmp_delta * spin(g, old_to) * (old_to_d ** params['phi'] - (old_to_d - 1) ** params['phi'])
    tmp_delta = 0.0

    # new_from neighbours
    neigs = g.neighbors(new_from)
    neigs_d = g.degree(neigs)
    for i, neig in enumerate(neigs):
        if neig not in (old_from, old_to, new_from, new_to):
            tmp_delta += (neigs_d[i] ** params['phi']) * spin(g, neig)
    delta += tmp_delta * spin(g, new_from) * (new_from_d ** params['phi'] - (new_from_d + 1) ** params['phi'])
    tmp_delta = 0.0

    # new_to neighbours
    neigs = g.neighbors(new_to)
    neigs_d = g.degree(neigs)
    for i, neig in enumerate(neigs):
        if neig not in (old_from, old_to, new_from, new_to):
            tmp_delta += (neigs_d[i] ** params['phi']) * spin(g, neig)
    delta += tmp_delta * spin(g, new_to) * (new_to_d ** params['phi'] - (new_to_d + 1) ** params['phi'])
    tmp_delta = 0.0

    # cross terms
    if (old_from, new_from) in edge_list or (new_from, old_from) in edge_list:
        delta += spin(g, old_from) * spin(g, new_from) * ((old_from_d ** params['phi']) * (new_from_d ** params['phi'])
                                                          - ((old_from_d - 1) ** params['phi']) * ((new_from_d + 1) ** params['phi']))

    if (old_from, new_to) in edge_list or (new_to, old_from) in edge_list:
        delta += spin(g, old_from) * spin(g, new_to) * ((old_from_d ** params['phi']) * (new_to_d ** params['phi'])
                                                          - ((old_from_d - 1) ** params['phi']) * ((new_to_d + 1) ** params['phi']))

    if (old_to, new_from) in edge_list or (new_from, old_to) in edge_list:
        delta += spin(g, old_to) * spin(g, new_from) * ((old_to_d ** params['phi']) * (new_from_d ** params['phi'])
                                                          - ((old_to_d - 1) ** params['phi']) * ((new_from_d + 1) ** params['phi']))

    if (old_to, new_to) in edge_list or (new_to, old_to) in edge_list:
        delta += spin(g, old_to) * spin(g, new_to) * ((old_to_d ** params['phi']) * (new_to_d ** params['phi'])
                                                        - ((old_to_d - 1) ** params['phi']) * ((new_to_d + 1) ** params['phi']))

    return delta


def energy_change_edge_gamma(g, old_from, old_to, new_from, new_to, params=None):
    old_from_d, old_to_d, new_from_d, new_to_d = g.degree([old_from, old_to, new_from, new_to])
    delta = old_from_d ** params['gamma'] - (old_from_d - 1) ** params['gamma'] \
            + old_to_d ** params['gamma'] - (old_to_d - 1) ** params['gamma'] \
            + new_from_d ** params['gamma'] - (new_from_d + 1) ** params['gamma'] \
            + new_to_d ** params['gamma'] - (new_to_d + 1) ** params['gamma']
    return delta


def new_vertex(n, exclude, params=None):
    new_v = np.int_(np.floor(rand() * n))
    if new_v in exclude:
        return new_vertex(n, exclude)
    return new_v


def draw_new_edge(n, edge_list, v_from, v_to, params=None):
    new_from = new_vertex(n, (v_from, v_to))
    new_to = new_vertex(n, (v_from, v_to, new_from))
    if (new_from, new_to) in edge_list or (new_to, new_from) in edge_list:
        return draw_new_edge(n, edge_list, v_from, v_to)
    return new_from, new_to


def one_step(g, n, m, params=None):
    v_index = np.int_(np.floor(rand() * n))
    e_index = np.int_(np.floor(rand() * m))
    changed = False

    # spin switching
    neigs = g.neighbors(v_index)
    delta = energy_change_spin(g, -2 * g.vs(v_index)["spin"][0][0], v_index, neigs, params=params)
    if delta <= 0.0 or rand() < np.exp(- params['beta'] * delta):
        g.vs(v_index)["spin"][0][0] *= -1
        g.vs(v_index)["color"] = GREEN if g.vs(v_index)["spin"][0][0] == 1 else RED
        changed = True

    # edge rewiring
    edge_list = g.get_edgelist()
    v_from, v_to = edge_list[e_index]
    new_from, new_to = draw_new_edge(n, edge_list, v_from, v_to, params=params)

    delta = energy_change_edge_gamma(g, v_from, v_to, new_from, new_to, params=params)
    delta += energy_change_edge_alpha(g, edge_list, v_from, v_to, new_from, new_to, params=params)
    if delta <= 0.0 or rand() < np.exp(- params['beta'] * delta):
        g.delete_edges([(v_from, v_to)])
        g.add_edges([(new_from, new_to)])
        g.es()['width'] = edge_width
        changed = True
    return g, changed


def count_incompatible(g, params=None):
    number = 0
    for v_from, v_to in g.get_edgelist():
        if spin(g, v_from) != spin(g, v_to):
            number += 1
    return number


def compute_energy(g, params=None):
    energy = 0.0
    for v_from, v_to in g.get_edgelist():
        energy += - ((g.degree(v_from) * g.degree(v_to)) ** params['phi']) * spin(g, v_from) * spin(g, v_to)

    for k in g.degree():
        energy += - (k ** params['gamma'])
    return energy


def main_loop(g, n, m, params=None):
    mag, mag_abs, energy, inc, c_num, largest_c, k_max = [], [], [], [], [], [], []
    k_dist = OrderedDict.fromkeys([k for k in xrange(0, n)], 0.0)  # min degree is 0, maximal is n-1

    for i in xrange(params['thermal']):
        g = one_step(g, n, m)
        if i % max(n, m) == 0:
            print('Thermalization: {}'.format(i))

    for i in xrange(params['sample_time']):
        g = one_step(g, n, m)
        if i % max(n, m) == 0:
            print('Sampling: {}'.format(i))
            mag.append(np.sum(g.vs()["spin"]))
            mag_abs.append(np.abs(mag[-1]))
            energy.append(compute_energy(g))
            inc.append(count_incompatible(g))
            c_num.append(len(g.clusters()))
            largest_c.append(len(max(g.clusters(), key=lambda clust: len(clust))))
            k_max.append(max(g.degree()))
            for x, _, y in g.degree_distribution().bins():
                k_dist[x] += y

    return g, np.array(mag, dtype=float) / n, np.array(mag_abs, dtype=float) / n, \
           np.array(energy, dtype=float) / (n * m), np.array(inc, dtype=float) / m, \
           np.array(c_num, dtype=float) / n, np.array(largest_c, dtype=float) / n, np.array(k_max, dtype=float) / n, \
           k_dist


def plot_degree_dist(params=None):
    graph = initialize_graph(params['n'], params['m'])

    start = time.time()
    graph, mag, mag_abs, energy, inc, c_num, largest_c, k_max, k_dist = main_loop(graph, params['n'], params['m'])
    end = time.time()

    print('Magnetization: {} +- {}'.format(np.abs(np.mean(mag)), np.std(mag)))
    print('Magnetization abs: {} +- {}'.format(np.mean(mag_abs), np.std(mag_abs)))
    print('Energy: {} +- {}'.format(np.mean(energy), np.std(energy)))
    print('Incompatible links: {} +- {}'.format(np.mean(inc), np.std(inc)))
    print('Largest component: {} +- {}'.format(np.mean(largest_c), np.std(largest_c)))
    print('Number of components: {} +- {}'.format(np.mean(c_num), np.std(c_num)))
    print('Maximal degree: {} +- {}'.format(np.mean(k_max), np.std(k_max)))

    with open('data/k_dist_N{}_M{}_T{}_GA{}_PH{}.csv'.format(params['n'], params['m'], params['temp'], params['gamma'], params['phi']), 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['k'] + [k for k in k_dist.keys()])
        writer.writerow(['y'] + [num for num in k_dist.values()])

    graph.write_pickle('data/graph_N{}_M{}_T{}_GA{}_PH{}.ig'.format(params['n'], params['m'], params['temp'], params['gamma'], params['phi']))

    print("Main loop finished in {} min.".format(round((end - start)/60.0), 1))


def plot_thermalization(g, n, m, times, params=None):
    mag, mag_abs, energy, inc, c_num, largest_c, k_max, corr, time_step = [], [], [], [], [], [], [], [], []

    for i in xrange(times):
        g = one_step(g, n, m)
        if i % 100 == 0:
            print('Step: {}'.format(i))
            mag.append(np.sum(g.vs()["spin"]))
            mag_abs.append(np.abs(mag[-1]))
            energy.append(compute_energy(g))
            inc.append(count_incompatible(g))
            c_num.append(len(g.clusters()))
            largest_c.append(len(max(g.clusters(), key=lambda clust: len(clust))))
            k_max.append(max(g.degree()))
            corr.append(g.assortativity_degree(directed=False))
            time_step.append(i)

    mag_abs, = plt.plot(time_step, np.array(mag_abs, dtype=float) / n, label='m_abs')
    inc, = plt.plot(time_step, np.array(inc, dtype=float) / m, label='inc')
    c_num, = plt.plot(time_step, np.array(c_num, dtype=float) / n, label='N_c')
    largest_c, = plt.plot(time_step, np.array(largest_c, dtype=float) / n, label='S')
    k_max, = plt.plot(time_step, np.array(k_max, dtype=float) / n, label='k_max')
    corr, = plt.plot(time_step, corr, label='corr')

    handles = [mag_abs, inc, c_num, largest_c, k_max, corr]
    plt.legend(handles=handles, fontsize=7)
    plt.xlabel('time step')
    plt.ylabel('order param')
    plt.savefig('plots/therm_N{}_M{}_GA{}_paramsPHI{}_T{}.png'.format(n, m, params['gamma'], params['phi'], params['temp']), format='png')
    plt.clf()

    plt.plot(time_step, np.array(energy, dtype=float) / (n * m))
    plt.xlabel('time step')
    plt.ylabel('energy')
    plt.savefig('plots/energy_therm_N{}_M{}_GA{}_PHI{}_T{}.png'.format(n, m, params['gamma'], params['phi'], params['temp']), format='png')
    plt.clf()

    for i in xrange(len(g.vs())):
        if g.vs(i)['spin'][0][0] == 1:
            g.vs(i)['color'] = '#26A57C'
        else:
            g.vs(i)['color'] = '#DA3A49'
    g.es["width"] = 1.5
    g.vs["size"] = 16
    lay = ig.Graph.layout_kamada_kawai(g)
    ig.plot(g, target='plots/graph_N{}_M{}_GA{}_PHI{}_T{}.png'.format(n, m, params['gamma'], params['phi'], params['temp']))
    ig.plot(g, layout=lay, target='plots/graph_kk_N{}_M{}_GA{}_PHI{}_T{}.png'.format(n, m, params['gamma'], params['phi'], params['temp']))


def plot_therm_multi(params=None):
    times = 100000
    params['gamma'] = 1.0
    params['phi'] = 0.0

    for t in np.linspace(0.1, 60, 8):
        for f in np.linspace(0, 2, 8):
            params['temp']= t
            params['beta'] = 1.0 / t
            params['phi'] = f
            graph = initialize_graph(params['n'], params['m'])
            plot_thermalization(graph, params['n'], params['m'], times)

    params['gamma'] = 1.0
    params['phi'] = 0.0

    for t in np.linspace(0.1, 60, 8):
        for g in np.linspace(0, 3, 8):
            params['temp']= t
            params['beta'] = 1.0 / t
            params['gamma'] = g
            graph = initialize_graph(params['n'], params['m'])
            plot_thermalization(graph, params['n'], params['m'], times)


if __name__ == '__main__':
    # plot_degree_dist()
    # plot_therm_multi()

    graph = initialize_graph(params['n'], params['m'])
    graph, mag, mag_abs, energy, inc, c_num, largest_c, k_max, k_dist = main_loop(graph, params['n'], params['m'])
    # lay = ig.Graph.layout(g)
    graph.write_pickle('data/graph_N{}_M{}_T{}_GA{}_PH{}.ig'.format(params['n'], params['m'],params['temp'], params['gamma'], params['phi']))
    lay = ig.Graph.layout_kamada_kawai(graph)
    ig.plot(graph, layout=lay)
