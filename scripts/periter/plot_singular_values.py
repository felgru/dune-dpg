#!/usr/bin/python
# -*- coding: utf8 -*-
import fileinput
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

if len(sys.argv) < 2:
    outname = 'singular_values.pdf'
else:
    outname = sys.argv[1]
henyeyGreenstein = re.compile(
            r'Henyey Greenstein kernel with gamma = ([0-9\.e\-\+]+)')
singularValues = list()
labels = list()
for line in fileinput.input(sys.argv[2:]):
    if fileinput.isfirstline():
        singularValues.append(list())
        m = henyeyGreenstein.search(line)
        if m:
            gamma = float(m.group(1))
            label = '$\gamma = {}$'.format(gamma)
        else:
            label = line
        labels.append(label)
    else:
        singularValues[-1].append(float(line))
numSingularValues = len(singularValues[0])
colors = ['#0054AF', # RWTH Blau
          '#57AB27', # RWTH Grün
          '#E30066', # RWTH Magenta
          '#BDCD00', # RWTH Maigrün
         ]
if len(singularValues) == 1:
    plt.plot(singularValues[0], linewidth=2.0, color=colors[0])
else:
    for c, (sv, l) in enumerate(zip(singularValues, labels)):
        plt.plot(sv, label=l, color=colors[c % len(colors)])
    plt.legend(loc='right')
plt.xlim(0, numSingularValues-1)
plt.yscale('log')
plt.savefig(outname)
