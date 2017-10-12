#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import fileinput
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import sys

if len(sys.argv) < 2:
    outname = 'conditions.pdf'
else:
    outname = sys.argv[1]
levels = list()
conditions = list()
line_pattern = re.compile(r'([0-9]*): ([0-9]*\.?[0-9]*e?[+-]?[0-9]?[0-9]?)')
for line in fileinput.input(sys.argv[2:]):
    m = line_pattern.match(line)
    levels.append(int(m.group(1)))
    conditions.append(float(m.group(2)))
plt.plot(levels, conditions)
plt.plot(levels, map(lambda l: 2**(2*l), levels), '.r')
plt.yscale('log')
plt.ylabel('condition of trial-to-test matrix')
plt.xlabel('level $\ell$, $h = 2^{-\ell}$')
plt.savefig(outname)
