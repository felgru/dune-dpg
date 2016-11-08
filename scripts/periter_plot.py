#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def readData(datafile):
    parametersPattern = re.compile(
        r'^Periter with rho = ([0-9]*\.?[0-9]*)'
        r', CT = ([0-9]*\.?[0-9]*)'
        r', kappa1 = ([0-9]*\.?[0-9]*)'
        r', kappa2 = ([0-9]*\.?[0-9]*)'
        r', kappa3 = ([0-9]*\.?[0-9]*)'
        , re.MULTILINE)
    dataPattern = re.compile(
        r'^Iteration ([0-9]+)\.([0-9]+): [^0-9]*([0-9\.e\-\+]+)'
        r', grid level: ([0-9]+), number of DOFs: ([0-9]+)'
        r', applying the kernel took ([0-9]+)us, (.*)$',
        re.MULTILINE)
    svdPattern = re.compile(
        r'SVD approximation with rank ([0-9]+)')
    waveletSVDPattern = re.compile(
        r'Wavelet SVD approximation with rank ([0-9]+) and level ([0-9]+)')
    waveletCompressionPattern = re.compile(
        r'MatrixCompression approximation with level ([0-9]+)')
    iterationIndices = list()
    gridResolutions = list()
    dofs = list()
    aposterioriErrors = list()
    kernelTimings = list()
    ranks = list()
    with open(datafile,"r") as errors:
        errors = errors.read()
        parametersMatch = parametersPattern.search(errors)
        parameters = { 'rho':    parametersMatch.group(1)
                     , 'CT':     parametersMatch.group(2)
                     , 'kappa1': parametersMatch.group(3)
                     , 'kappa2': parametersMatch.group(4)
                     , 'kappa3': parametersMatch.group(5)
                     }
        for (n, nRefinement, aPostErr, gridLevel, numDOFs, time, rest) \
                in dataPattern.findall(errors):
            iterationIndices.append(n+'.'+nRefinement)
            n = int(n)
            nRefinement = int(nRefinement)
            aPostErr = float(aPostErr)
            gridLevel = int(gridLevel)
            numDOFs = int(numDOFs)
            time = int(time) / 1000000.;
            gridResolutions.append(np.exp2(-gridLevel))
            dofs.append(numDOFs)
            aposterioriErrors.append(aPostErr)
            kernelTimings.append(time)
            m = svdPattern.match(rest)
            if m:
                ranks.append(m.group(1))
            else:
                m = waveletSVDPattern.match(rest)
                if m:
                    ranks.append(m.group(1))
                else:
                    m = waveletCompressionPattern.match(rest)
                    if m:
                        ranks.append('-')
    return { 'parameters': parameters
           , 'datapoints': len(iterationIndices)
           , 'iterationIndices': iterationIndices
           , 'gridResolutions': gridResolutions
           , 'dofs': dofs
           , 'aposterioriErrors': aposterioriErrors
           , 'kernelTimings': kernelTimings
           , 'ranks': ranks
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
    print((r'convergence table for $\rho = {p[rho]}$'
           r', $C_T = {p[CT]}$, $\kappa_1 = {p[kappa1]}$'
           r', $\kappa_2 = {p[kappa2]}$, $\kappa_3 = {p[kappa3]}$'
           '\n'
          ).format(p=data['parameters']))
    print(r'\begin{tabular}{l|lllll}')
    print('iteration & $h$ & \#DOFs & aposteriori error & duration of kernel application / s & rank of kernel approximation \\\\\n\hline')
    for i in range(data['datapoints']):
        d = { 'iterationIndex': data['iterationIndices'][i]
            , 'gridResolution': data['gridResolutions'][i]
            , 'dofs': data['dofs'][i]
            , 'aposterioriError': data['aposterioriErrors'][i]
            , 'kernelTiming': data['kernelTimings'][i]
            , 'rank': data['ranks'][i]
            }
        print((r'{d[iterationIndex]} & {d[gridResolution]} & {d[dofs]}'
               r'& {d[aposterioriError]} & {d[kernelTiming]} & {d[rank]}\\'
              ).format(d=d))
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
