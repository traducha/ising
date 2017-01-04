import sys
import csv
import matplotlib.pyplot as plt

# f_name = "{}.csv".format(sys.argv[1])
x_max = 20000000 # int(sys.argv[2])

for iii in []:#[1,2]:
		f_name = "res/energy_vs_B_{}.csv".format(iii)
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i])
		    std[i] = float(std[i])

		plt.errorbar(x[1:x_max], value[1:x_max], yerr=std[1:x_max], fmt='o')
		plt.title(f_name)
		plt.xlabel(x[0])
		plt.ylabel(value[0])
		plt.savefig(f_name.replace("csv", "png"), format="png")
		plt.clf()

		f_name = "res/mag_vs_B_{}.csv".format(iii)
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i])
		    std[i] = float(std[i])


if 0:
		f_name = "res/mag_vs_B_1.csv"
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i])
		    std[i] = float(std[i])

		plt.scatter(x[1:x_max], value[1:x_max], color='r')
		plt.title(f_name)
		plt.xlabel(x[0])
		plt.ylabel(value[0])

		f_name = "res/mag_vs_B_2.csv"
		with open(f_name, 'rb') as file:
		    x, value, std = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(x)):
		    x[i] = float(x[i])
		    value[i] = float(value[i])
		    std[i] = float(std[i])

		plt.scatter(x[1:x_max], value[1:x_max])
		plt.title("1 - red, 2 - blue")
		plt.xlabel(x[0])
		plt.ylabel(value[0])
		plt.savefig("algorithms_comparsion_mag.png", format="png")
		plt.clf()


