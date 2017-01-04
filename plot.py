import sys
import csv
import matplotlib.pyplot as plt

# f_name = "{}.csv".format(sys.argv[1])
x_max = 20000000 # int(sys.argv[2])

for iii in [2]:#[1,2]:
	for j in range(2):
		f_name = "res/energy{}_{}.csv".format(iii, j)
		with open(f_name, 'rb') as file:
		    index, value = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(index)):
		    index[i] = int(index[i])
		    value[i] = float(value[i])

		plt.plot(index[1:x_max], value[1:x_max])
		plt.title(f_name)
		plt.xlabel(index[0])
		plt.ylabel(value[0])
		plt.savefig(f_name.replace("csv", "png"), format="png")
		plt.clf()

		f_name = "res/mag{}_{}.csv".format(iii, j)
		with open(f_name, 'rb') as file:
		    index, value = csv.reader(file, delimiter=';', quotechar='|')

		for i in xrange(1, len(index)):
		    index[i] = int(index[i])
		    value[i] = float(value[i])

		plt.plot(index[1:x_max], value[1:x_max])
		plt.title(f_name)
		plt.xlabel(index[0])
		plt.ylabel(value[0])
		plt.savefig(f_name.replace("csv", "png"), format="png")
		plt.clf()


