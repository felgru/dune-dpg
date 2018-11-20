#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import argparse
import fileinput
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import sys

# to fix cut-off x label
mpl.rcParams.update({'figure.autolayout': True})

aparser = argparse.ArgumentParser(description='Plot condition of trial-to-test matrix')
aparser.add_argument('-o', dest='outfile', default='conditions.pdf')
aparser.add_argument('--presentation', dest='presentation',
                     action='store_true', default=False,
                     help='generate plot for presentations')
args, unk = aparser.parse_known_args()
levels = list()
conditions = list()
line_pattern = re.compile(r'([0-9]*): ([0-9]*\.?[0-9]*e?[+-]?[0-9]?[0-9]?)')
for line in fileinput.input(unk):
    m = line_pattern.match(line)
    levels.append(int(m.group(1)))
    conditions.append(float(m.group(2)))
if args.presentation:
    mpl.rc('font',**{'family':'serif','size':20})
    mpl.rc('text', usetex=True)
    plt.plot(levels, conditions, linewidth=2.0, color='#0054AF')
    plt.plot(levels, map(lambda l: 2**(2*l), levels), '.r', linewidth=2.0)
else:
    plt.plot(levels, conditions)
    plt.plot(levels, map(lambda l: 2**(2*l), levels), '.r')
plt.yscale('log')
plt.ylabel('condition of trial-to-test matrix')
plt.xlabel('level $\ell$, $h = 2^{-\ell}$')
plt.savefig(args.outfile)
