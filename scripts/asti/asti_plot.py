#!/usr/bin/python3
# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function)
import argparse
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from output_parser import readData

def plot_convergence(data,
         outputfile='asti_error.pdf',
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

    rhoN = [ (float(data['params']['rho']))**k for k in np.arange(len(data['globalAccIteratesDiff']))]
    errIdealIteration = []
    for n in range(len(rhoN)):
        t = np.asarray(data['eta'])[0:n+1][::-1]
        errIdealIteration.append(np.sum(rhoN[0:n+1]*t))

    iterationIndices = data['iterationIndices']

    line1 = ax1.plot(iterationIndices, data['accKernel'],
                      label='$k_n$: err kernel approx')

    line1_ = ax1.plot(iterationIndices, data['aPost'],
                     label='$t_n$: err transport solves (a posteriori estimation)')

    if not simple_plot:
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
        plt.setp(line1___, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[3])
        plt.setp(line1____, linewidth=2.0,
                 marker='o', markersize=4.0,
                 color=colorPalette[4])
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
         outputfile='asti_error.pdf',
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
         outputfile='asti_error.pdf',
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
         outputfile='asti_error.pdf',
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
         outputfile='asti_error.pdf',
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
         outputfile='asti_error.pdf',
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
            '#cc9900', '#ffd24d',  # yellow
            ],
         simple_plot=False):
    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(0,0))

    rhoN = [ (float(data['params']['rho']))**k for k in np.arange(len(data['globalAccIteratesDiff']))]
    errIdealIteration = []
    for n in range(len(rhoN)):
        t = np.asarray(data['eta'])[0:n+1][::-1]
        errIdealIteration.append(np.sum(rhoN[0:n+1]*t))

    iterationIndices = data['iterationIndices']
    innerIterationStats = data['innerIterationsStats']
    minNumInnerIterations = []
    maxNumInnerIterations = []
    avgNumInnerIterations = []
    for oi in iterationIndices:
        minNum = sys.maxsize
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
        minNum = sys.maxsize
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

def _plot_Dofs_per_direction(data,
         ax,
         ylabel='#DoFs / direction',
         ylim=None,
         yscale='log',
         color='#0063cc'):
    ax.set_ylabel(ylabel, color=color)
    ax.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    innerIterationsStats = data['innerIterationsStats']
    num_Dofs_per_iteration = [ num_Dofs_per_direction(innerIterationsStats[oi])
                                for oi in iterationIndices ]

    pointPos = []
    pointVal = []
    violinPos = []
    violinVal = []
    for i, num_Dofs in zip(iterationIndices, num_Dofs_per_iteration):
        minNum = sys.maxsize
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
        pointPlot = ax.plot(pointPos, pointVal, 'o',
                            color=color)
    if violinPos:
        violinPlot = ax.violinplot(violinVal,
                                   positions=violinPos,
                                   showmeans=True,
                                   showmedians=False)
        plt.setp(violinPlot['bodies'], color=color)
        plt.setp(violinPlot['cmeans'], color=color)
        plt.setp(violinPlot['cmins'], color=color)
        plt.setp(violinPlot['cmaxes'], color=color)
        plt.setp(violinPlot['cbars'], color=color)
        #plt.setp(violinPlot['cmedians'], color=color)

    ax.set_yscale(yscale)
    if ylim != None:
        ax.set_ylim(ylim)

def plot_Dofs_per_direction(data,
         outputfile='asti_dofs.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel='#DoFs / direction',
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='log',
         colorPalette=[
            '#0063cc', '#80bdff',  # blue
            '#33cc33', '#99e699',  # green
            '#cc0000', '#ff5c33',  # red
            '#b800e6', '#e580ff',  # purple
            '#cc9900', '#ffd24d'  # yellow
            ]):
    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)

    _plot_Dofs_per_direction(data,
         ax,
         ylabel=ylabel,
         ylim=ylim,
         yscale=yscale,
         color=colorPalette[0])

    ax.set_xscale(xscale)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    plt.savefig(outputfile)

    plt.clf()

def _plot_num_directions(data,
         ax,
         ylabel='#directions',
         ylim=None,
         yscale='linear',
         color='#cc0000'):
    ax.set_ylabel(ylabel, color=color)
    ax.ticklabel_format(style='sci', scilimits=(0,0))

    iterationIndices = data['iterationIndices']
    innerIterationsStats = data['innerIterationsStats']
    num_Dofs_per_iteration = [ num_Dofs_per_direction(innerIterationsStats[oi])
                                for oi in iterationIndices ]

    directions = ax.plot(iterationIndices, [len(nd) for nd in num_Dofs_per_iteration])
    # plot in RWTH red
    plt.setp(directions, linewidth=2.0,
             marker='x', markersize=4.0,
             color=color)

    ax.set_yscale(yscale)
    if ylim != None:
        ax.set_ylim(ylim)

def plot_num_directions(data,
         outputfile='asti_directions.pdf',
         title=None,
         xlabel='Outer Iteration',
         ylabel='#directions',
         xlim=None,
         ylim=None,
         xscale='linear',
         yscale='linear',
         colorPalette=[
            '#0063cc', '#80bdff',  # blue
            '#33cc33', '#99e699',  # green
            '#cc0000', '#ff5c33',  # red
            '#b800e6', '#e580ff',  # purple
            '#cc9900', '#ffd24d'  # yellow
            ]):
    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)

    _plot_num_directions(data,
         ax,
         ylabel=ylabel,
         ylim=ylim,
         yscale=yscale,
         color=colorPalette[4])

    ax.set_xscale(xscale)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    plt.savefig(outputfile)

    plt.clf()

def plot_Dofs_and_directions_vs_iteration(data,
         outputfile='asti_dofs.pdf',
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
            ]):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    if title != None:
        plt.title(title)
    ax1.set_xlabel(xlabel)

    _plot_Dofs_per_direction(data,
         ax1,
         ylabel=ylabel[0],
         ylim=ylim[0] if ylim is not None else None,
         yscale=yscale[0] if yscale is not None else None,
         color=colorPalette[0])

    _plot_num_directions(data,
         ax2,
         ylabel=ylabel[1],
         ylim=ylim[1] if ylim is not None else None,
         yscale=yscale[1] if yscale is not None else None,
         color=colorPalette[4])

    ax1.set_xscale(xscale)
    ax1.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    if xlim != None:
        plt.xlim(xlim)
    plt.savefig(outputfile)

    plt.clf()

def plot_a_posteriori_err_VS_dofs(data,
         outputfile='asti_a_posteriori_VS_dofs.pdf',
         title=None,
         xlabel='#DoFs',
         ylabel='a posteriori error',
         xlim=None,
         ylim=None,
         xscale='log',
         yscale='log',
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

    # plot in RWTH blue
    plt.plot('dofs', 'globalAccAposteriori', 'o-', data=data,
             color=colorPalette[0])

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    ax.grid()
    plt.savefig(outputfile)

    plt.clf()

def parse_ylim(ylim_string):
    if not ylim_string:
        return None
    return [float(y) for y in ylim_string.split(',')]

def parse_ylims(ylim_string):
    if not ylim_string:
        return None
    ylims = parse_ylim(ylim_string)
    return (ylims[0:2], ylims[2:])


aparser = argparse.ArgumentParser(
        description='Generate convergence plot and table for ASTI')
#aparser.add_argument('--preamble', dest='print_preamble',
#                     action='store_true', default=False,
#                     help='print Latex preamble for the convergence table')
aparser.add_argument('--paper', dest='simple_plot',
                     action='store_true', default=False,
                     help='generate simpler plot for our paper')
aparser.add_argument('--combined-dofs-and-directions-plot',
                     dest='combined_plot',
                     action='store_true', default=False,
                     help='combine plots of #DoFs and #directions')
aparser.add_argument('--conv-ylim', dest='conv_ylim', action='store',
                     help='limits of the y-axis of the convergence plot')
aparser.add_argument('--dofs-ylim', dest='dofs_ylim', action='store',
                     help='limits of the y-axes of the dofs plot')
aparser.add_argument('--dirs-ylim', dest='dirs_ylim', action='store',
                     help='limits of the y-axes of the directions plot')
aparser.add_argument('--apost-lim', dest='apost_lim', action='store',
                     help='limits for the axes of the a posteriori / dofs plot')
aparser.add_argument('infile', action='store')
aparser.add_argument('prefixOutputFile', action='store',
                     help='prefix of the name of the plot files')
args = aparser.parse_args()
args.conv_ylim = parse_ylim(args.conv_ylim)
args.dofs_ylim = parse_ylim(args.dofs_ylim)
args.dirs_ylim = parse_ylim(args.dirs_ylim)
args.apost_lim = parse_ylims(args.apost_lim)
if args.apost_lim:
    args.apost_xlim = args.apost_lim[0]
    args.apost_ylim = args.apost_lim[1]

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
     simple_plot=args.simple_plot,
     ylim=args.conv_ylim
    )
plot_directions(data,
     outputfile=args.prefixOutputFile+"-directions.pdf",
    )

plot_kernel_acc_VS_time(data,
     outputfile=args.prefixOutputFile+"-kernel-acc-VS-time.pdf",
    )

plot_kernel_matrix_info(data,
     outputfile=args.prefixOutputFile+"-kernel-matrix-info.pdf",
    )

if(data['params']['kernelApproxType'] == 'SVD'):
    plot_svd(data,
     outputfile=args.prefixOutputFile+"-svd.pdf",
    )

plot_inner_iterations(data,
     outputfile=args.prefixOutputFile+"-inner-iterations.pdf",
     simple_plot=args.simple_plot
    )

if args.combined_plot:
    plot_Dofs_and_directions_vs_iteration(data,
         outputfile=args.prefixOutputFile+"-num-dofs.pdf",
         ylim=(args.dofs_ylim, args.dirs_ylim)
        )
else:
    plot_Dofs_per_direction(data,
         outputfile=args.prefixOutputFile+"-num-dofs.pdf",
         ylim=args.dofs_ylim
        )
    plot_num_directions(data,
         outputfile=args.prefixOutputFile+"-num-directions.pdf",
         ylim=args.dirs_ylim
        )

plot_a_posteriori_err_VS_dofs(data,
     outputfile=args.prefixOutputFile+"-a-posteriori-VS-dofs.pdf",
     xlim=args.apost_xlim,
     ylim=args.apost_ylim,
    )
