#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def readData(datafile):
    dataPattern = re.compile(
        '^Iteration ([0-9]+)\.([0-9]+): [^0-9]*([0-9\.e\-\+]+)'
        ', grid level: ([0-9]+)',
        re.MULTILINE)
    gridResolutions = list()
    aposterioriErrors = list()
    with open(datafile,"r") as errors:
        errors = errors.read()
        for (n, nRefinement, aPostErr, gridLevel) \
                in dataPattern.findall(errors):
            n = int(n)
            nRefinement = int(nRefinement)
            aPostErr = float(aPostErr)
            gridLevel = int(gridLevel)
            gridResolutions.append(np.exp2(-gridLevel))
            aposterioriErrors.append(aPostErr)
    return (gridResolutions, aposterioriErrors)

def plot(gridResolutions,
         errors,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='mesh size $H$',
         ylabel='a posteriori error estimators',
         xlim=None,
         ylim=None,
         xscale='log',
         yscale='log'):
    fig = plt.figure()
    if title != None:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    line = plt.plot(gridResolutions, errors)
    # plot in RWTH blue
    plt.setp(line, linewidth=2.0,
             marker='o', markersize=3.0,
             color='#0054AF')

    plt.xscale(xscale)
    plt.yscale(yscale)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()



if len(sys.argv) != 3:
    print('Usage: sys.argv[0] infile outfile')
    sys.exit(1)

(gridResolutions, errors) = readData(sys.argv[1])

mpl.rc('text', usetex=True)

plot(gridResolutions,
     errors,
     outputfile=sys.argv[2],
     # title='a posteriori errors of Periter',
    )
