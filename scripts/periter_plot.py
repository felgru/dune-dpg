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
        r'^Periter with a priori accuracy( up to)? ([0-9]*\.?[0-9]*)'
        r' for the kernel'
        r', rho = ([0-9]*\.?[0-9]*)'
        r', CT = ([0-9]*\.?[0-9]*)'
        r', kappa1 = ([0-9]*\.?[0-9]*)'
        r', kappa2 = ([0-9]*\.?[0-9]*)'
        r', kappa3 = ([0-9]*\.?[0-9]*)'
        , re.MULTILINE)
    dataPattern = re.compile(
        r'^Error at end of Iteration ([0-9]+): ([0-9]+\.?[0-9]*)'
        r', using ([0-9]+) DoFs'
        r', accuracy was ([0-9]+\.?[0-9]*)'
        r', eta was ([0-9]+\.?[0-9]*)'
        r', applying the kernel took ([0-9]+)us, (.*)$',
        re.MULTILINE)
    svdPattern = re.compile(
        r'SVD approximation with rank ([0-9]+)')
    haarWaveletSVDPattern = re.compile(
        r'Wavelet SVD approximation with rank ([0-9]+) and level ([0-9]+)')
    haarWaveletCompressionPattern = re.compile(
        r'MatrixCompression approximation with level ([0-9]+)')
    alpertWaveletSVDPattern = re.compile(
        r'Wavelet SVD approximation with Alpert wavelets of order ([0-9]+)'
        r', rank ([0-9]+) and level ([0-9]+)')
    iterationIndices = list()
    dofs = list()
    targetAccuracies = list()
    etas = list()
    aposterioriErrors = list()
    kernelTimings = list()
    ranks = list()
    with open(datafile,"r") as errors:
        errors = errors.read()
        parametersMatch = parametersPattern.search(errors)
        parameters = { 'adaptiveInS': parametersMatch.group(1) != ''
                     , 'kernelAccuracy': parametersMatch.group(2)
                     , 'rho':    parametersMatch.group(3)
                     , 'CT':     parametersMatch.group(4)
                     , 'kappa1': parametersMatch.group(5)
                     , 'kappa2': parametersMatch.group(6)
                     , 'kappa3': parametersMatch.group(7)
                     }
        for (n, aPostErr, numDOFs, targetAccuracy, eta, time, rest) \
                in dataPattern.findall(errors):
            iterationIndices.append(int(n))
            dofs.append(int(numDOFs))
            targetAccuracies.append(float(targetAccuracy))
            etas.append(float(eta))
            aposterioriErrors.append(float(aPostErr))
            kernelTimings.append(int(time) / 1000000.)
            m = svdPattern.match(rest)
            if m:
                ranks.append(m.group(1))
            else:
                m = haarWaveletSVDPattern.match(rest)
                if m:
                    ranks.append(m.group(1))
                else:
                    m = alpertWaveletSVDPattern.match(rest)
                    if m:
                        ranks.append(m.group(2))
                    else:
                        m = haarWaveletCompressionPattern.match(rest)
                        if m:
                            ranks.append('-')
    return { 'parameters': parameters
           , 'datapoints': len(iterationIndices)
           , 'iterationIndices': iterationIndices
           , 'dofs': dofs
           , 'targetAccuracies': targetAccuracies
           , 'etas': etas
           , 'aposterioriErrors': aposterioriErrors
           , 'kernelTimings': kernelTimings
           , 'ranks': ranks
           }

def plot(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='outer iteration',
         ylabel=('a posteriori error estimator', '# of DoFs'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best'):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel[0])
    ax2.set_ylabel(ylabel[1])
    ax1.ticklabel_format(style='sci', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    line1 = ax1.plot(iterationIndices, data['aposterioriErrors'],
                     label='a posteriori error')
    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=3.0,
             color='#0054AF')
    line1_ = ax1.plot(iterationIndices, data['targetAccuracies'],
                      label='estimate for '
                            '$\\rho^nC_{\mathcal{T}}\|f\|_{V\'}+2\eta_n$')
    line1__ = ax1.plot(iterationIndices, data['etas'], label='$\eta$')

    line2 = ax2.plot(iterationIndices, data['dofs'], label='# of DoFs')
    # plot in RWTH purple
    plt.setp(line2, linewidth=2.0,
             marker='x', markersize=3.0,
             color='#612158')

    ax1.set_xscale(xscale)
    ax2.set_xscale(xscale)
    ax1.set_yscale(yscale)
    ax2.set_yscale(yscale)
    if legendlocation != None:
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        plt.legend(lines1 + lines2, labels1 + labels2,
                   loc=legendlocation, shadow=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()

def print_table(data):
    if data['parameters']['adaptiveInS']:
        up_to = 'up to'
    else:
        up_to = ''
    print((r'convergence table for $\rho = {p[rho]}$'
           r', $C_T = {p[CT]}$, $\kappa_1 = {p[kappa1]}$'
           r', $\kappa_2 = {p[kappa2]}$, $\kappa_3 = {p[kappa3]}$'
           r' with accuracy {up_to} {p[kernelAccuracy]} in the'
           r' kernel approximation\n'
          ).format(p=data['parameters'], up_to=up_to))
    print(r'\begin{tabular}{r|rrrl}')
    print(r'& \multicolumn{2}{c}{kernel approximation} & & \\')
    print('iteration & duration / s & rank & \#DOFs & aposteriori error \\\\\n')
    for row in range(len(data['iterationIndices'])):
        print(r'\hline')
        print(r'{n} '.format(n=data['iterationIndices'][row]))
        print(r'& {t} & {r} '
                .format(t=data['kernelTimings'][row], r=data['ranks'][row]))
        print(r'& {dofs} & {err} \\'
               .format(dofs=data['dofs'][row],
                        err=data['aposterioriErrors'][row]))
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

#mpl.rc('text', usetex=True)

plot(data,
     outputfile=args.outfile,
     # title='a posteriori errors of Periter',
    )
