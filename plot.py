import sys
import csv
import matplotlib.pyplot as plt

f_name = "magnetization.csv"

with open(f_name, 'rb') as file:
    index, value = csv.reader(file, delimiter=';', quotechar='|')

for i in xrange(1, len(index)):
    index[i] = int(index[i])
    value[i] = int(value[i])

plt.plot(index[1:], value[1:])
plt.title(f_name)
plt.xlabel(index[0])
plt.ylabel(value[0])
plt.show()

