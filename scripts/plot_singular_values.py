#!/usr/bin/python
import fileinput
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    outname = 'singular_values.pdf'
else:
    outname = sys.argv[1]
singularValues = list()
for line in fileinput.input(sys.argv[2:]):
    singularValues.append(float(line))
plt.plot(singularValues)
plt.xlim(0, len(singularValues)-1)
plt.yscale('log')
plt.savefig(outname)
