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
    # parametersPattern = re.compile(
    #     r'^Periter with( up to)? ([0-9]+) directions, rho = ([0-9]*\.?[0-9]*)'
    #     r', CT = ([0-9]*\.?[0-9]*)'
    #     r', kappa1 = ([0-9]*\.?[0-9]*)'
    #     r', kappa2 = ([0-9]*\.?[0-9]*)'
    #     r', kappa3 = ([0-9]*\.?[0-9]*)'
    #     , re.MULTILINE)
    parametersPattern = re.compile(
        r'^PERITER algorithm\n'
        r'=================\n'
        r'Prescribed final accuracy: ([0-9]*\.?[0-9]*)\n'
        r'Henyey Greenstein kernel with gamma = ([0-9]*\.?[0-9]*)\n'
        r'Wavelet order: ([0-9]*\.?[0-9]*)\n'
        r'Kernel approximation with: (.*)\n'
        r'Maximum wavelet level: ([0-9]*\.?[0-9]*)\n'
        r'Maximum number of directions: ([0-9]*\.?[0-9]*)\n'
        r'Periter parameters:\n'
        r'rho = ([0-9]*\.?[0-9]*)\n'
        r'rhobar = ([0-9]*\.?[0-9]*)\n'
        r'kappa1 = ([0-9]*\.?[0-9]*)\n'
        r'kappa2 = ([0-9]*\.?[0-9]*)\n'
        r'kappa3 = ([0-9]*\.?[0-9]*)\n'
        r'CT = ([0-9]*\.?[0-9]*)\n'
        , re.MULTILINE)
    iterationIndicesPattern = re.compile(r'Iteration n=([0-9]*\.?[0-9]*)\n')
    etaPattern = re.compile(r'eta_n = rhobar\^{-n}: ([0-9]*\.?[0-9]*)\n')
    wltLevelPattern = re.compile(r'Current wavelet level: ([0-9]*\.?[0-9]*)\n')
    numSPattern = re.compile(r'Number of directions: ([0-9]*\.?[0-9]*)\n')
    svdRankPattern = re.compile(r'SVD rank: ([0-9]*\.?[0-9]*)\n')
    matrixTHpattern = re.compile(
        r'Kernel matrix is of size ([0-9]*\.?[0-9]*)x([0-9]*\.?[0-9]*).'
        r' It has ([0-9]*\.?[0-9]*) elements'
        r' of which ([0-9]*\.?[0-9]*) are zero.\n'
        , re.MULTILINE)
    timeEvalKernelPattern = re.compile(r'Computing time: ([0-9]*\.?[0-9]*)us')
    aPostPattern = re.compile(r'Total a posteriori error: ([0-9]*\.?[0-9]*)\n')
    accKernelPattern = re.compile(r'Accuracy kernel: ([0-9]*\.?[0-9]*)\n')
    globalAccApostPattern = re.compile(
        r'Global accuracy \(a posteriori\): ([0-9]*\.?[0-9]*)\n')
    globalAccAprioriPattern = re.compile(
        r'Global accuracy \(a priori\): ([0-9]*\.?[0-9]*)')
    dofsPattern = re.compile(r'Total number of DoFs: ([0-9]*\.?[0-9]*)\n')
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
        parameters = { 'eps': parametersMatch.group(1)
                     , 'gamma':   parametersMatch.group(2)
                     , 'wltOrder':    parametersMatch.group(3)
                     , 'kernelApproxType':    parametersMatch.group(4)
                     , 'maxWltLevel':     parametersMatch.group(5)
                     , 'maxNumS': parametersMatch.group(6)
                     , 'rho': parametersMatch.group(7)
                     , 'rhobar': parametersMatch.group(8)
                     , 'kappa1': parametersMatch.group(9)
                     , 'kappa2': parametersMatch.group(10)
                     , 'kappa3': parametersMatch.group(11)
                     , 'CT': parametersMatch.group(12)
                     }
        iterationIndices = iterationIndicesPattern.findall(errors)
        eta = etaPattern.findall(errors)
        wltLevel = wltLevelPattern.findall(errors)
        numS = numSPattern.findall(errors)
        svdRank = svdRankPattern.findall(errors)
        matrixTH = matrixTHpattern.findall(errors)
        timeEvalKernel = timeEvalKernelPattern.findall(errors)
        aPost = aPostPattern.findall(errors)
        accKernel = accKernelPattern.findall(errors)
        globalAccApost = globalAccApostPattern.findall(errors)
        globalAccApriori = globalAccAprioriPattern.findall(errors)
        dofs = dofsPattern.findall(errors)

    return { 'params': parameters
           , 'iterationIndices': iterationIndices
           , 'eta': eta
           , 'wltLevel': wltLevel
           , 'numS': numS
           , 'svdRank': svdRank
           , 'matrixTH': matrixTH
           , 'timeEvalKernel': timeEvalKernel
           , 'aPost': aPost
           , 'accKernel': accKernel
           , 'globalAccApost': globalAccApost
           , 'globalAccApriori': globalAccApriori
           , 'dofs': dofs
           }

def plot_convergence(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel=('Error','# DoFs'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=['#0054AF','#612158','#33cc33','#cc3300','#cc9900']):
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
    line1 = ax1.plot(iterationIndices, data['aPost'],
                     label='Error Transport (a posteriori estimation)')

    line1_ = ax1.plot(iterationIndices, data['accKernel'],
                      label='Error Kernel approx')

    line1_ = ax1.plot(iterationIndices, data['globalAccApost'],
                      label='Global Error (Transport+Kernel)')

    line1__ = ax1.plot(iterationIndices, data['eta'], label='$\eta_n$')

    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[1])
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[2])
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[3])

    line2 = ax2.plot(iterationIndices, data['dofs'], label='# of DoFs')
    # plot in RWTH purple
    plt.setp(line2, linewidth=2.0,
             marker='x', markersize=4.0,
             color=colorPalette[4])

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

def plot_directions(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel=(''),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='linear',
         legendlocation='best',
         colorPalette=['#0054AF','#612158','#33cc33','#cc3300','#cc9900']):
    fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    # ax2.set_ylabel(ylabel[1])
    ax1.ticklabel_format(style='sci', scilimits=(0,0))
    # ax2.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    line1 = ax1.plot(iterationIndices, data['numS'],
                     label='# directions')

    line1_ = ax1.plot(iterationIndices, data['wltLevel'],
                      label='Wavelet level')

    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])
    plt.setp(line1_, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[1])


    # line2 = ax2.plot(iterationIndices, data['numS'], label='# directions')
    # # plot in RWTH purple
    # plt.setp(line2, linewidth=2.0,
    #          marker='x', markersize=3.0,
    #          color='#612158')

    ax1.set_xscale(xscale)
    # ax2.set_xscale(xscale)
    ax1.set_yscale(yscale)
    # ax2.set_yscale(yscale)
    if legendlocation != None:
        lines1, labels1 = ax1.get_legend_handles_labels()
        # lines2, labels2 = ax2.get_legend_handles_labels()
        plt.legend(lines1, labels1,
                   loc=legendlocation, shadow=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()

def plot_kernel_acc_VS_time(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='outer iteration',
         ylabel=('Accuracy kernel','Computing time kernel eval'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=['#0054AF','#612158','#33cc33','#cc3300','#cc9900']):
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
    line1 = ax1.plot(iterationIndices, data['accKernel'],
                      label='Error Kernel approx')
    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])

    line2 = ax2.plot(iterationIndices, data['timeEvalKernel'],
                      label='Computing time for kernel evaluation')
    # plot in RWTH purple
    plt.setp(line2, linewidth=2.0,
             marker='x', markersize=4.0,
             color=colorPalette[1])

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

def plot_kernel_matrix_info(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='outer iteration',
         ylabel=('Accuracy kernel','Computing time kernel eval (in $\mu$s)'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=['#0054AF','#612158','#33cc33','#cc3300','#cc9900']):
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
    if(data['params']['kernelApproxType'] == 'SVD'):
        line1 = ax1.plot(iterationIndices, data['svdRank'],
                      label='SVD rank')
    else:
        if(data['params']['kernelApproxType'] == 'Matrix compression'):
            totalentriesKernelMatrix \
                = np.asarray([c[2] for c in data['matrixTH']], dtype=float)
            zerosKernelMatrix = np.asarray([c[3] for c in data['matrixTH']], dtype=float)
            ratio = zerosKernelMatrix/totalentriesKernelMatrix
            line1 = ax1.plot(iterationIndices, ratio,
                      label='# zeros entries / # entries in kernel matrix')
        else:
            print('Function plot_kernel_matrix_info: unsupported kernelApproxType')

    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])

    line2 = ax2.plot(iterationIndices, data['timeEvalKernel'],
                    label='Computing time kernel evaluation (in $\mu$s)')
    # plot in RWTH purple
    plt.setp(line2, linewidth=2.0,
             marker='x', markersize=4.0,
             color=colorPalette[1])

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

# TODO: Adapt to new version
def print_table(data):
    if data['parameters']['adaptiveInS']:
        up_to = 'up to'
    else:
        up_to = ''
    print((r'convergence table for $\rho = {p[rho]}$'
           r', $C_T = {p[CT]}$, $\kappa_1 = {p[kappa1]}$'
           r', $\kappa_2 = {p[kappa2]}$, $\kappa_3 = {p[kappa3]}$'
           r' with {up_to} {p[numS]} directions'
           '\n'
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

# TODO: Adapt to new version
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
aparser.add_argument('prefixOutputFile', action='store',
                     help='prefix of the name of the plot files')
args = aparser.parse_args()

data = readData(args.infile)

# # TODO: Adapt to new version
# if args.print_preamble:
#     print_preamble()
#     print_table(data)
#     print(r'\end{document}')
# else:
#     print_table(data)

#mpl.rc('text', usetex=True)

plot_convergence(data,
     outputfile=args.prefixOutputFile+"-conv.pdf",
     # title='a posteriori errors of Periter',
    )
plot_directions(data,
     outputfile=args.prefixOutputFile+"-directions.pdf",
     # title='a posteriori errors of Periter',
    )
plot_kernel_acc_VS_time(data,
     outputfile=args.prefixOutputFile+"-kernel-acc-VS-time.pdf",
     # title='a posteriori errors of Periter',
    )

plot_kernel_matrix_info(data,
     outputfile=args.prefixOutputFile+"-kernel-matrix-info.pdf",
     # title='a posteriori errors of Periter',
    )
