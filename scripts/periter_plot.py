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
        ', grid level: ([0-9]+)'
        ', applying the kernel took ([0-9]+)us, (.*)$',
        re.MULTILINE)
    iterationIndices = list()
    gridResolutions = list()
    aposterioriErrors = list()
    kernelTiming = list()
    with open(datafile,"r") as errors:
        errors = errors.read()
        for (n, nRefinement, aPostErr, gridLevel, time, rest) \
                in dataPattern.findall(errors):
            iterationIndices.append(n+'.'+nRefinement)
            n = int(n)
            nRefinement = int(nRefinement)
            aPostErr = float(aPostErr)
            gridLevel = int(gridLevel)
            time = int(time) / 1000000.;
            gridResolutions.append(np.exp2(-gridLevel))
            aposterioriErrors.append(aPostErr)
            kernelTiming.append(time)
    return { 'datapoints': len(iterationIndices)
           , 'iterationIndices': iterationIndices
           , 'gridResolutions': gridResolutions
           , 'aposterioriErrors': aposterioriErrors
           , 'kernelTiming': kernelTiming
           }

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

def print_table(data):
    print(r'\begin{tabular}{l|lll}')
    print('iteration & $h$ & aposteriori error & duration of kernel application / s \\\\\n\hline')
    for i in range(data['datapoints']):
        d = { 'iterationIndex': data['iterationIndices'][i]
            , 'gridResolution': data['gridResolutions'][i]
            , 'aposterioriError': data['aposterioriErrors'][i]
            , 'kernelTiming': data['kernelTiming'][i]
            }
        print(r'{d[iterationIndex]} & {d[gridResolution]} & {d[aposterioriError]} & {d[kernelTiming]} \\'.format(d=d))
    print(r'\end{tabular}')



if len(sys.argv) != 3:
    print('Usage: sys.argv[0] infile outfile')
    sys.exit(1)

data = readData(sys.argv[1])

print_table(data)

mpl.rc('text', usetex=True)

plot(data['gridResolutions'],
     data['aposterioriErrors'],
     outputfile=sys.argv[2],
     # title='a posteriori errors of Periter',
    )
