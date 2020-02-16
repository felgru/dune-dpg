#!/usr/bin/python3
# -*- coding: utf8 -*-
import argparse
import fileinput
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

aparser = argparse.ArgumentParser(description='Plot profiling data')
aparser.add_argument('-o', dest='outfile', default='profile.pdf')
args, unk = aparser.parse_known_args()

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':16})
mpl.rc('text', usetex=True)

line_pattern = re.compile(r'([0-9]*)\s*([0-9]*)\s*([0-9]*)\s*([0-9]*)')
h = list()
buffered = list()
unbuffered = list()
umfpack = list()
for line in fileinput.input(unk):
    match = line_pattern.match(line)
    h.append(1./float(match.group(1)))
    buffered.append(float(match.group(2))/1000000.)
    unbuffered.append(float(match.group(3))/1000000.)
    umfpack.append(float(match.group(4))/1000000.)

plt.xlabel('mesh size $H$')
plt.ylabel('runtime in s')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xscale('log')
plt.yscale('log')
plt.plot(h, buffered, label='buffered TestspaceCoefficientMatrix',
         marker='o')
plt.plot(h, unbuffered, label='unbuffered TestspaceCoefficientMatrix',
         marker='x')
plt.plot(h, umfpack, label='solving the system with UMFPACK')
plt.legend(loc='lower left', shadow=True, fontsize='small')
plt.savefig(args.outfile)
