#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function)
import argparse
import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict

def readData(datafile):
    parametersPattern = re.compile(
        r'PERITER algorithm\n'
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
    singularValuesPattern = re.compile(
        r'Singular values of kernel matrix:\n'
        r'(([0-9]+\.?[0-9]*e?-?[0-9]*)\n)*'
        , re.MULTILINE)
    iterationIndicesPattern = re.compile(r'Iteration n=([0-9]+)\n')
    etaPattern = re.compile(r'eta_n = rhobar\^{-n}: ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
    wltLevelPattern = re.compile(r'Current wavelet level: ([0-9]+)\n')
    numSPattern = re.compile(r'Number of directions: ([0-9]+)\n')
    svdRankPattern = re.compile(r'SVD rank: ([0-9]+)\n')
    matrixTHpattern = re.compile(
        r'Kernel matrix is of size ([0-9]+)x([0-9]+).'
        r' It has ([0-9]+) elements'
        r' of which ([0-9]+) are zero.\n'
        , re.MULTILINE)
    timeEvalKernelPattern = re.compile(r'Computing time: ([0-9]*\.?[0-9]*)us')
    aPostPattern = re.compile(r'Error transport solves \(a posteriori estimation\): ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
    accKernelPattern = re.compile(r'Accuracy kernel: ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
    globalAccIterationApostPattern = re.compile(
        r'Error bound \|\|bar u_n -T\^{-1}K bar u_{n-1}\|\| \(a posteriori\): ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
    globalAccIteratesDiffPattern = re.compile(r'Error bound \|\|u_n - bar u_n\|\| \(a posteriori\): ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
    globalAccAprioriPattern = re.compile(
        r'A priori bound global accuracy \|\|u - bar u_n\|\|: ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)')
    globalAccAposterioriPattern = re.compile(
        r'A posteriori bound global accuracy \|\|u - bar u_n\|\|: ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)')
    dofsPattern = re.compile(r'Total number of DoFs: ([0-9]+)\n')
    innerIterationsPattern = re.compile(
            r'Iteration ([0-9]+)\.([0-9]+) for direction [0-9]+:\n'
            r'  - A posteriori estimation of \|\| \(u,trace u\) - '
                  r'\(u_fem,theta\) \|\| = [0-9]*\.?[0-9]*e?[+-]?[0-9]+?\n'
            r'  - Grid level: ([0-9]+)\n'
            r'  - Number of DOFs: ([0-9]+)\n\n'
            r'a posteriori error for current direction: [0-9]*\.?[0-9]*e?[+-]?[0-9]+? '
                  r'\((enough|not enough)'
            , re.MULTILINE)
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
        singularValues = []
        if(parameters['kernelApproxType']=='SVD'):
            svPat = re.compile(r'([0-9]+\.?[0-9]*e?-?[0-9]*)\n', re.MULTILINE)
            singularValues = svPat.findall(singularValuesPattern.search(errors).group())
        iterationIndices = iterationIndicesPattern.findall(errors)
        etas = etaPattern.findall(errors)
        wltLevel = wltLevelPattern.findall(errors)
        numS = numSPattern.findall(errors)
        svdRank = svdRankPattern.findall(errors)
        matrixTH = matrixTHpattern.findall(errors)
        timeEvalKernel = timeEvalKernelPattern.findall(errors)
        aPost = aPostPattern.findall(errors)
        accKernel = accKernelPattern.findall(errors)
        globalAccIterationApost = globalAccIterationApostPattern.findall(errors)
        globalAccIteratesDiff = globalAccIteratesDiffPattern.findall(errors)
        globalAccApriori = globalAccAprioriPattern.findall(errors)
        globalAccAposteriori = globalAccAposterioriPattern.findall(errors)
        dofs = dofsPattern.findall(errors)
        innerIterationsStats = defaultdict(list)
        for m in innerIterationsPattern.finditer(errors):
            innerIterationsStats[m.group(1)].append(
                    { 'numIterations': int(m.group(2))
                    , 'maxLevel': int(m.group(3))
                    , 'numDOFs': int(m.group(4))
                    })

    return { 'params': parameters
           , 'singularValues': singularValues
           , 'iterationIndices': iterationIndices
           , 'eta': etas
           , 'wltLevel': wltLevel
           , 'numS': numS
           , 'svdRank': svdRank
           , 'matrixTH': matrixTH
           , 'timeEvalKernel': timeEvalKernel
           , 'aPost': aPost
           , 'accKernel': accKernel
           , 'globalAccIterationApost': globalAccIterationApost
           , 'globalAccIteratesDiff' : globalAccIteratesDiff
           , 'globalAccApriori': globalAccApriori
           , 'globalAccAposteriori': globalAccAposteriori
           , 'dofs': dofs
           , 'innerIterationsStats': innerIterationsStats
           }

def plot_convergence(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel='Error',
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='lower left',
         colorPalette=[
            '#0063cc', '#80bdff',  # blue
            '#33cc33', '#99e699',  # green
            '#cc0000', '#ff5c33',  # red
            '#b800e6', '#e580ff',  # purple
            '#cc9900', '#ffd24d'  # yellow
            ],
         simple_plot=False):
    fig, ax1 = plt.subplots()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.ticklabel_format(style='sci', scilimits=(0,0))

    rhoN = [ (float(data['params']['rho']))**k for k in np.arange(len(data['globalAccIterationApost']))]
    errIdealIteration = []
    for n in range(len(rhoN)):
        t = ((np.asarray(map(float, data['eta'])))[0:n+1])[::-1]
        errIdealIteration.append(np.sum(rhoN[0:n+1]*t))

    iterationIndices = data['iterationIndices']

    line1 = ax1.plot(iterationIndices, data['accKernel'],
                      label='$k_n$: err kernel approx')

    line1_ = ax1.plot(iterationIndices, data['aPost'],
                     label='$t_n$: err transport solves (a posteriori estimation)')

    if not simple_plot:
        line1__ = ax1.plot(iterationIndices
                        , data['globalAccIterationApost']
                        , label=r'$e_n = t_n+C_T k_n$ ($||\bar u_n -T^{-1}K'
                                r'\bar u_{n-1}||\leq e_n)$')

        line1___ = ax1.plot(iterationIndices, data['eta'],
                            label=r'$\eta_n (e_n\leq\eta_n)$')

    line1____ = ax1.plot(iterationIndices, data['globalAccIteratesDiff'],
        label=r'a posteriori bound for $||u_n - \bar u_n||$')

    lineAposteriori = ax1.plot(iterationIndices, data['globalAccAposteriori'],
        label=r'a posteriori bound for $||u - \bar u_n||$')

    if not simple_plot:
        line1_____ = ax1.plot(iterationIndices, errIdealIteration,
            label=r'$\sum_{j=0}^{n} \rho^j \eta_{n-j}$ '
                  r'($\sum_{j=0}^{n} \rho^j e_{n-j} '
                  r'\leq \sum_{j=0}^{n} \rho^j \eta_{n-j}$)')

        line1______ = ax1.plot(iterationIndices,
                (1.+np.pi*np.pi/6.)*np.asarray(rhoN),
                label=r'$(1+\pi^2/6)\rho^n$ '
                      r'($\sum_{j=0}^{n} \rho^j \eta_{n-j} '
                      r'\leq (1+\pi^2/6)\rho^n$)')

    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])
    plt.setp(line1_, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[1])
    if not simple_plot:
        plt.setp(line1__, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[2])
        plt.setp(line1___, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[3])
    plt.setp(line1____, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[4])
    if not simple_plot:
        plt.setp(line1_____, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[5])
        plt.setp(line1______, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[6])
    plt.setp(lineAposteriori, linewidth=2.0,
             marker='x', markersize=4.0,
             color=colorPalette[7])

    ax1.set_xscale(xscale)
    ax1.set_yscale(yscale)
    # Shrink current axis by 20%
    if not simple_plot:
        box1 = ax1.get_position()
        ax1.set_position([box1.x0, box1.y0,
            box1.width, box1.height * 0.6])
    if legendlocation != None:
        lines1, labels1 = ax1.get_legend_handles_labels()
        if simple_plot:
            plt.legend(lines1, labels1,
                       loc=legendlocation, shadow=True,
                       ncol=1, fancybox=True,fontsize=12)
        else:
            plt.legend(lines1, labels1,
                       loc=legendlocation, shadow=True,
                       bbox_to_anchor=(0.5, 1.9),
                       ncol=1, fancybox=True, fontsize=12)
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
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
         ylabel='',
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='linear',
         legendlocation='best',
         colorPalette=[
         '#0063cc', '#80bdff',  # blue
         '#33cc33', '#99e699',  # green
         '#cc0000', '#ff5c33',  # red
         '#b800e6', '#e580ff',  # purple
         '#cc9900', '#ffd24d'  # yellow
         ]):
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
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()

def plot_svd(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel=(''),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=[
         '#0063cc', '#80bdff',  # blue
         '#33cc33', '#99e699',  # green
         '#cc0000', '#ff5c33',  # red
         '#b800e6', '#e580ff',  # purple
         '#cc9900', '#ffd24d'  # yellow
         ]):
    fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    # ax2.set_ylabel(ylabel[1])
    ax1.ticklabel_format(style='sci', scilimits=(0,0))
    # ax2.ticklabel_format(style='sci', scilimits=(0,0))

    line1 = ax1.plot(data['singularValues'],
                     label='Singular values of kernel matrix')

    # plot in RWTH blue
    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])

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
         ylabel=('Error Kernel approx',
            'Computing time kernel eval (in $\mu$s)'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=[
         '#0063cc', '#80bdff',  # blue
         '#33cc33', '#99e699',  # green
         '#cc0000', '#ff5c33',  # red
         '#b800e6', '#e580ff',  # purple
         '#cc9900', '#ffd24d'  # yellow
         ]):
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
                      label='Computing time for kernel evaluation (in $\mu$s)')
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
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
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
         ylabel=(('SVD rank','# zeros entries / # entries in kernel matrix'),
            'Computing time kernel eval (in $\mu$s)'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         legendlocation='best',
         colorPalette=[
         '#0063cc', '#80bdff',  # blue
         '#33cc33', '#99e699',  # green
         '#cc0000', '#ff5c33',  # red
         '#b800e6', '#e580ff',  # purple
         '#cc9900', '#ffd24d'  # yellow
         ]):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel[1])
    ax1.ticklabel_format(style='sci', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    if(data['params']['kernelApproxType'] == 'SVD'):
        ax1.set_ylabel(ylabel[0][0])
        line1 = ax1.plot(iterationIndices, data['svdRank'],
                      label='SVD rank')
    else:
        if(data['params']['kernelApproxType'] == 'Matrix compression'):
            ax1.set_ylabel(ylabel[0][1])
            totalentriesKernelMatrix \
                = np.asarray([c[2] for c in data['matrixTH']], dtype=float)
            zerosKernelMatrix = np.asarray([c[3] for c in data['matrixTH']], dtype=float)
            ratio = zerosKernelMatrix/totalentriesKernelMatrix
            line1 = ax1.plot(iterationIndices, ratio,
                      label='# zeros entries / # entries in kernel matrix')
        else:
            print('Function plot_kernel_matrix_info: unsupported kernelApproxType')

    plt.setp(line1, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])

    line2 = ax2.plot(iterationIndices, data['timeEvalKernel'],
                    label='Computing time kernel evaluation (in $\mu$s)')
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
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()

def plot_inner_iterations(data,
         outputfile='periter_error.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel='inner iterations',
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='linear',
         legendlocation='upper center',
         colorPalette=[
            '#0063cc', '#80bdff',  # blue
            '#33cc33', '#99e699',  # green
            '#cc0000', '#ff5c33',  # red
            '#b800e6', '#e580ff',  # purple
            '#cc9900', '#ffd24d'  # yellow
            ],
         simple_plot=False):
    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(0,0))

    rhoN = [ (float(data['params']['rho']))**k for k in np.arange(len(data['globalAccIterationApost']))]
    errIdealIteration = []
    for n in range(len(rhoN)):
        t = ((np.asarray(map(float, data['eta'])))[0:n+1])[::-1]
        errIdealIteration.append(np.sum(rhoN[0:n+1]*t))

    iterationIndices = data['iterationIndices']
    innerIterationStats = data['innerIterationsStats']
    minNumInnerIterations = []
    maxNumInnerIterations = []
    avgNumInnerIterations = []
    for oi in iterationIndices:
        minNum = sys.maxint
        maxNum = 0
        sum = 0
        for direction in innerIterationStats[oi]:
            md = direction['numIterations']
            minNum = min(minNum, md)
            maxNum = max(maxNum, md)
            sum += md
        minNumInnerIterations.append(minNum)
        maxNumInnerIterations.append(maxNum)
        avgNumInnerIterations.append(sum / len(innerIterationStats[oi]))

    minLevel = []
    maxLevel = []
    avgLevel = []
    for oi in iterationIndices:
        minNum = sys.maxint
        maxNum = 0
        sum = 0
        for direction in innerIterationStats[oi]:
            md = direction['maxLevel']
            minNum = min(minNum, md)
            maxNum = max(maxNum, md)
            sum += md
        minLevel.append(minNum)
        maxLevel.append(maxNum)
        avgLevel.append(sum / len(innerIterationStats[oi]))

    minLine = ax.plot(iterationIndices, minNumInnerIterations,
                      label='min. number of inner iterations')

    avgLine = ax.plot(iterationIndices, avgNumInnerIterations,
                      label='avg. number of inner iterations')

    maxLine = ax.plot(iterationIndices, maxNumInnerIterations,
                      label='max. number of inner iterations')

    minLevelLine = ax.plot(iterationIndices, minLevel,
                      label='min. grid level')

    avgLevelLine = ax.plot(iterationIndices, avgLevel,
                      label='avg. grid level')

    maxLevelLine = ax.plot(iterationIndices, maxLevel,
                      label='max. grid level')

    # plot in RWTH blue
    plt.setp(minLine, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[0])

    plt.setp(avgLine, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[1])

    plt.setp(maxLine, linewidth=2.0,
             marker='o', markersize=4.0,
             color=colorPalette[2])

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    # Shrink current axis by 20%
    if legendlocation != None:
        lines1, labels1 = ax.get_legend_handles_labels()
        if simple_plot:
            plt.legend(lines1, labels1,
                       loc=legendlocation, shadow=True,
                       ncol=1, fancybox=True,fontsize=12)
        else:
            plt.legend(lines1, labels1,
                       loc=legendlocation, shadow=True,
                       bbox_to_anchor=(0.5, 1.9),
                       ncol=1, fancybox=True, fontsize=12)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()

def num_Dofs_per_direction(innerIterationStats):
    num_Dofs = [ it['numDOFs'] for it in innerIterationStats ]
    return num_Dofs

def plot_Dofs_per_direction(data,
         outputfile='periter_dofs.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel=('#DoFs / direction', '#directions'),
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale=('log', 'linear'),
         colorPalette=[
            '#0063cc', '#80bdff',  # blue
            '#33cc33', '#99e699',  # green
            '#cc0000', '#ff5c33',  # red
            '#b800e6', '#e580ff',  # purple
            '#cc9900', '#ffd24d'  # yellow
            ],
         simple_plot=False):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel[0], color=colorPalette[0])
    ax2.set_ylabel(ylabel[1], color=colorPalette[4])
    ax1.ticklabel_format(style='sci', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    innerIterationsStats = data['innerIterationsStats']
    num_Dofs_per_iteration = [ num_Dofs_per_direction(innerIterationsStats[oi])
                                for oi in iterationIndices ]

    pointPos = []
    pointVal = []
    violinPos = []
    violinVal = []
    for i, num_Dofs in enumerate(num_Dofs_per_iteration):
        minNum = sys.maxint
        maxNum = 0
        for num in num_Dofs:
            minNum = min(minNum, num)
            maxNum = max(maxNum, num)
        if minNum == maxNum:
            pointPos.append(i)
            pointVal.append(maxNum)
        else:
            violinPos.append(i)
            violinVal.append(num_Dofs)
    # plot in RWTH blue
    if pointPos:
        pointPlot = ax1.plot(pointPos, pointVal, 'o',
                             color=colorPalette[0])
    if violinPos:
        violinPlot = ax1.violinplot(violinVal,
                                    positions=violinPos,
                                    showmeans=True,
                                    showmedians=False)
        plt.setp(violinPlot['bodies'], color=colorPalette[0])
        plt.setp(violinPlot['cmeans'], color=colorPalette[0])
        plt.setp(violinPlot['cmins'], color=colorPalette[0])
        plt.setp(violinPlot['cmaxes'], color=colorPalette[0])
        plt.setp(violinPlot['cbars'], color=colorPalette[0])
        #plt.setp(violinPlot['cmedians'], color=colorPalette[0])

    directions = ax2.plot(iterationIndices, map(len, num_Dofs_per_iteration))
    # plot in RWTH red
    plt.setp(directions, linewidth=2.0,
             marker='x', markersize=4.0,
             color=colorPalette[4])

    ax1.set_xscale(xscale)
    ax1.set_yscale(yscale[0])
    ax2.set_xscale(xscale)
    ax2.set_yscale(yscale[1])
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        ax1.set_ylim(ylim[0])
        ax2.set_ylim(ylim[1])
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

def parse_ylim(ylim_string):
    if not ylim_string:
        return None
    return map(float,ylim_string.split(','))

def parse_ylims(ylim_string):
    if not ylim_string:
        return None
    ylims = parse_ylim(ylim_string)
    return (ylims[0:2], ylims[2:])


aparser = argparse.ArgumentParser(
        description='Generate convergence plot and table for Periter')
#aparser.add_argument('--preamble', dest='print_preamble',
#                     action='store_true', default=False,
#                     help='print Latex preamble for the convergence table')
aparser.add_argument('--paper', dest='simple_plot',
                     action='store_true', default=False,
                     help='generate simpler plot for our paper')
aparser.add_argument('--conv-ylim', dest='conv_ylim', action='store',
                     help='limits of the y-axis of the convergence plot')
aparser.add_argument('--dofs-ylim', dest='dofs_ylim', action='store',
                     help='limits of the y-axes of the dofs plot')
aparser.add_argument('infile', action='store')
aparser.add_argument('prefixOutputFile', action='store',
                     help='prefix of the name of the plot files')
args = aparser.parse_args()
args.conv_ylim = parse_ylim(args.conv_ylim)
args.dofs_ylim = parse_ylims(args.dofs_ylim)

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
     simple_plot=args.simple_plot,
     ylim=args.conv_ylim
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

if(data['params']['kernelApproxType'] == 'SVD'):
    plot_svd(data,
     outputfile=args.prefixOutputFile+"-svd.pdf",
     # title='a posteriori errors of Periter',
    )

plot_inner_iterations(data,
     outputfile=args.prefixOutputFile+"-inner-iterations.pdf",
     simple_plot=args.simple_plot
    )

plot_Dofs_per_direction(data,
     outputfile=args.prefixOutputFile+"-num-dofs.pdf",
     simple_plot=args.simple_plot,
     ylim=args.dofs_ylim
    )
