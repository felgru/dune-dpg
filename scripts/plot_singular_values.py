#!/usr/bin/python
import fileinput
import matplotlib.pyplot as plt

singularValues = list()
for line in fileinput.input():
    singularValues.append(float(line))
plt.plot(singularValues)
plt.yscale('log')
plt.savefig('singular_values.pdf')
