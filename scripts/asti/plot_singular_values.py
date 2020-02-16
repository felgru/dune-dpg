#!/usr/bin/python3
# -*- coding: utf8 -*-
import argparse
import fileinput
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

aparser = argparse.ArgumentParser(description='Plot singular values')
aparser.add_argument('-o', dest='outfile', default='singular_values.pdf')
aparser.add_argument('--presentation', dest='presentation',
                     action='store_true', default=False,
                     help='generate plot for presentations')
args, unk = aparser.parse_known_args()

henyeyGreenstein = re.compile(
            r'Henyey Greenstein kernel with gamma = ([0-9\.e\-\+]+)')
singularValues = list()
labels = list()
for line in fileinput.input(unk):
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
if args.presentation:
    mpl.rc('font',**{'family':'serif','size':20})
    linewidth = 2.0
else:
    linewidth = 1.0
if len(singularValues) == 1:
    plt.plot(singularValues[0], linewidth=linewidth, color=colors[0])
else:
    for c, (sv, l) in enumerate(zip(singularValues, labels)):
        plt.plot(sv, label=l, linewidth=linewidth,
                 color=colors[c % len(colors)])
    plt.legend(loc='right')
plt.xlim(0, numSingularValues-1)
plt.yscale('log')
plt.savefig(args.outfile)
