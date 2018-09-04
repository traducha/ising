import sys
import csv
import re
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import pprint
import matplotlib as mpl
from mean_field_2a import mean_field_phi as phi_mf
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['font.family'] = 'serif'

N = 1000
M = 3000
G = '1.000000'
A = '0.700000'
temp_lim = [0, 150]
y_lim = [0, 2]
y_values = np.linspace(y_lim[0], y_lim[1], 81)[:-1]
y_var = 'A'  # 'A' or 'G'
aspect = temp_lim[1] / (y_lim[1] - y_lim[0])

quants = ['mag_abs', 'largest_degree', 'energy', 'incompatible', 'mag', 'largest_clust', 'clust_num', 'degree_corr']
norm = [1.0*N, 1.0*N, 1.0*N*M, 1.0*M, 1.0*N, 1.0*N, 1.0*N, 1.0]

os.chdir("../res_2D/phi")

phis = np.linspace(0.0, 2.0, 1000)
temp = np.linspace(0.1, 150, 1000)

trans_temp = []
first = True

for p in phis:
    _, degrees = phi_mf(float(N), float(M), p, temp)
    degrees = degrees * N
    t = 0.0
    for i, degree in enumerate(degrees):
        if degree < 50:
            t = temp[i]
            break
    # print p, t
    # print degrees
    trans_temp.append(t)

final_p, final_t = [], []
for i, t in enumerate(trans_temp):
    if 0.0 < t < 1000:
        final_t.append(t)
        final_p.append(phis[i])
print final_t
print final_p
print

fig = plt.figure(figsize=(7.85, 4))
for j, q in enumerate(quants):
    if q == 'mag_abs':
        ax1 = plt.subplot(121)
        ax = ax1
    elif q == 'largest_degree':
        ax2 = plt.subplot(122, sharey=ax1)
        ax = ax2
    else:
        continue

    value_matrix = []
    std_matrix = []

    for y in y_values:
        g = G
        a = A
        if y_var == 'A':
            a = str(round(y, 6)).ljust(8, '0')
        elif y_var == 'G':
            g = str(round(y, 6)).ljust(8, '0')
        else:
            raise

        f_name = "{}_vs_B_N{}_L{}_J1.000000_h0.000000_FI0.500000_GA{}_AL{}.csv".format(q, N, M, g, a)
        with open(f_name, 'rb') as file:
            x, value, std = csv.reader(file, delimiter=';', quotechar='|')

        for i in xrange(1, len(x)):
            x[i] = float(x[i])
            value[i] = float(value[i]) / norm[j]
            std[i] = float(std[i]) / norm[j]

        value_matrix.append(value[1:])
        std_matrix.append(std[1:])

    plt.xlabel(r'$T = 1/ \beta$', fontsize=18)
    ax.set_xlim([0.0, 150.0])
    print np.min(value_matrix), np.max(value_matrix)
    if q == 'mag_abs':
        ax.set_title(r'$|m|$', fontsize=18)
        ax.set_ylabel(r'$\phi$' if y_var == 'A' else r'$\gamma$', fontsize=22)
        im1 = ax.imshow(value_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none',
                       aspect=aspect/1.07,
                       vmin=0.0, vmax=1.0)
        ax.set_title(r'$|m|$', fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=13, length=4, width=1.1)
    elif q == 'largest_degree':
        ax.set_title(r'$k_{max}$', fontsize=18)
        im2 = ax.imshow(value_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none',
                       aspect=aspect/1.07,
                       vmin=0.0, vmax=1.0)
        plt.plot([0, 150], [1.23723724, 1.23723724], '--', color='black', linewidth=2)
        plt.plot(final_t, final_p, color='#FFFF33', linewidth=2)
        ax.tick_params(axis='both', which='major', labelsize=13, length=4, width=1.1)

#Create and remove the colorbar for the first subplot
cbar1 = fig.colorbar(im1, ax=ax1)
fig.delaxes(fig.axes[2])

#Create second colorbar
cbar2 = fig.colorbar(im2, ax=ax2)


plt.setp(ax2.get_yticklabels(), visible=False)
plt.tight_layout()
plt.subplots_adjust(wspace=-0.19)
# plt.show()
plt.savefig('/home/tomaszraducha/Pulpit/2D_phi.pdf', format='pdf')
plt.clf()

    # im = plt.imshow(std_matrix, cmap=None, origin='lower', extent=temp_lim + y_lim, interpolation='none', aspect=aspect)
    # plt.colorbar(im)
    # plt.xlabel(r'$T = 1/ \beta$', fontsize=16)
    # plt.ylabel(r'$\phi$' if y_var == 'A' else r'$\gamma$', fontsize=20)
    # plt.tight_layout()
    # plt.show()
    # # plt.savefig('/home/tomasz/Desktop/2D/gamma/g_STD_2D{}.pdf'.format(q), format='pdf')
    # plt.clf()
