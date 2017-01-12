#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import argparse
import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def readData(datafile):
    parametersPattern = re.compile(
        r'^Periter with ([0-9]+) directions, rho = ([0-9]*\.?[0-9]*)'
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
        parameters = { 'numS':   parametersMatch.group(1)
                     , 'rho':    parametersMatch.group(2)
                     , 'CT':     parametersMatch.group(3)
                     , 'kappa1': parametersMatch.group(4)
                     , 'kappa2': parametersMatch.group(5)
                     , 'kappa3': parametersMatch.group(6)
                     }
        for (n, nRefinement, aPostErr, gridLevel, numDOFs, time, rest) \
                in dataPattern.findall(errors):
            iterationIndices.append((n, nRefinement))
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
           r' with {p[numS]} directions'
           '\n'
          ).format(p=data['parameters']))
    data2 = list()
    for i in range(data['datapoints']):
        (it, it_inner) = data['iterationIndices'][i]
        d = { 'innerIteration': it_inner
            , 'gridResolution': data['gridResolutions'][i]
            , 'dofs': data['dofs'][i]
            , 'aposterioriError': data['aposterioriErrors'][i]
            }
        if data2 and (data2[-1]['n'] == it):
            data2[-1]['inner_data'].append(d)
        else:
            data2.append({ 'n': it
                         , 'kernelTiming': data['kernelTimings'][i]
                         , 'rank': data['ranks'][i]
                         , 'inner_data': [d]
                         })
    print(r'\begin{tabular}{rr|rrlrl}')
    print(r'& & \multicolumn{2}{c}{kernel approximation} & & & \\')
    print('\\multicolumn{2}{c|}{iteration} & duration / s & rank & $h$ & \#DOFs & aposteriori error \\\\\n')
    for row in data2:
        print(r'\hline')
        print(r'\multirow{{{s}}}{{*}}{{{n}}} '
                .format(n=row['n'], s=len(row['inner_data'])))
        print(r'& {} '.format(row['inner_data'][0]['innerIteration']))
        print(r'& \multirow{{{s}}}{{*}}{{{t}}} '
                .format(t=row['kernelTiming'], s=len(row['inner_data'])))
        print(r'& \multirow{{{s}}}{{*}}{{{r}}} '
                .format(r=row['rank'], s=len(row['inner_data'])))
        print((r'& {d[gridResolution]} '
               r'& {d[dofs]} & {d[aposterioriError]} \\'
              ).format(d=row['inner_data'][0]))
        for d in row['inner_data'][1:]:
            print((r'& {d[innerIteration]} & & & {d[gridResolution]} '
                   r'& {d[dofs]} & {d[aposterioriError]} \\'
                  ).format(d=d))
    print(r'\end{tabular}')

def print_preamble():
    print(r'\documentclass[11pt,a4paper]{article}' '\n'
          '%\n'
          '% For narrow margins\n'
          r'\usepackage{fullpage}' '\n'
          '\n'
          r'\usepackage[utf8]{inputenc}' '\n'
          '\n'
          '% packages from the American Mathematical Society (AMS)\n'
          r'\usepackage{amsmath, amsthm, amssymb}' '\n'
          '% More fancy functionality for theorems\n'
          r'\usepackage{thmtools, thm-restate}' '\n'
          '\n'
          r'% For \MoveEqLeft, \coloneqq, etc.' '\n'
          r'\usepackage{mathtools}' '\n'
          '\n'
          r'\usepackage{multirow}' '\n'
          '\n'
          r'\begin{document}')


aparser = argparse.ArgumentParser(
        description='Generate convergence plot and table for Periter')
aparser.add_argument('--preamble', dest='print_preamble',
                     action='store_true', default=False,
                     help='print Latex preamble for the convergence table')
aparser.add_argument('infile', action='store')
aparser.add_argument('outfile', action='store',
                     help='name of the convergence plot file')
args = aparser.parse_args()

data = readData(args.infile)

if args.print_preamble:
    print_preamble()
    print_table(data)
    print(r'\end{document}')
else:
    print_table(data)

mpl.rc('text', usetex=True)

plot(data['gridResolutions'],
     data['aposterioriErrors'],
     outputfile=args.outfile,
     # title='a posteriori errors of Periter',
    )
