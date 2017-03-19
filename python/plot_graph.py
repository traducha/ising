import igraph as ig
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from pprint import pprint
mpl.rcParams['font.family'] = 'serif'

N = 100
M = 300
T = 10.0
GAMMA = 2.0
PHI = 0.0


def poisson(k, lambda_):
    if isinstance(k, (list, tuple, np.ndarray)):
        return [1.0 * (lambda_ ** k_) * (np.exp(-lambda_)) / (math.factorial(k_)) for k_ in k]
    return 1.0 * (lambda_ ** k) * (np.exp(-lambda_)) / (math.factorial(k))


def plot_degree(temp=T, save=False):
    with open('k_dist_N{}_M{}_T{}_GA{}_PH{}.csv'.format(N, M, temp, GAMMA, PHI), 'rb') as _file:
        k, y = csv.reader(_file, delimiter=';', quotechar='|')

    for i in xrange(1, len(k)):
        k[i] = int(k[i])
        y[i] = float(y[i])
    print('Number of lonely nodes: {}'.format(y[1]))
    pprint(y)

    k = k[1:]
    y = np.array(y[1:]) / np.sum(y[1:])

    if 0:  # 1 to omit k=0
        k = k[1:]
        y = y[1:]

    plt.plot(k, poisson(k, 2.0 * M / N), color='#3078A8')
    plt.scatter(k, y, marker='o', facecolors='#ff5e00', edgecolors='#ff5e00')
    plt.title(r"Degree distribution for $N$={}, $M$={}, $T$={}, $\gamma$={}, $\phi$={}".format(N, M, temp, GAMMA, PHI))
    plt.xlabel(r"$k$")
    plt.ylabel(r"$P(k)$")
    # plt.xscale('log')
    # plt.yscale('log')
    if save:
        plt.savefig('k_dist_N{}_M{}_T{}_GA{}_PH{}.png'.format(N, M, temp, GAMMA, PHI), format="png")
    else:
        plt.show()
    plt.clf()


def plot_graph(temp=T, save=False):
    g = ig.Graph.Read_Pickle('graph_N{}_M{}_T{}_GA{}_PH{}.ig'.format(N, M, temp, GAMMA, PHI))
    for i in xrange(len(g.vs())):
        if g.vs(i)['spin'][0][0] == 1:
            g.vs(i)['color'] = '#26A57C'
        else:
            g.vs(i)['color'] = '#DA3A49'
    g.es["width"] = 1.5

    if save:
        # possible formats .png, .pdf, .ps
        ig.plot(g, target='graph_N{}_M{}_T{}_GA{}_PH{}.png'.format(N, M, temp, GAMMA, PHI))
    else:
        ig.plot(g)


if __name__ == '__main__':
    for t in (0.5, 10.0, 20.0, 30.0, 40.0):
        plot_degree(temp=t, save=True)
        plot_graph(temp=t, save=True)
