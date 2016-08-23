#!/usr/bin/python
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def insertData(data, ls, ks, la, ka, n, l2error, aposteriori):
    if ls not in data:
        data[ls] = dict()
    if ks not in data[ls]:
        data[ls][ks] = dict()
    if la not in data[ls][ks]:
        data[ls][ks][la] = dict()
    if ka not in data[ls][ks][la]:
        data[ls][ks][la][ka] = dict()
    if 'n' not in data[ls][ks][la][ka]:
        data[ls][ks][la][ka]['n']     = list()
        data[ls][ks][la][ka]['l2']    = list()
        data[ls][ks][la][ka]['apost'] = list()
    data[ls][ks][la][ka]['n'].append(n)
    data[ls][ks][la][ka]['l2'].append(l2error)
    data[ls][ks][la][ka]['apost'].append(aposteriori)

def readData(data, datafile):
    dataPattern = re.compile(
        '^([0-9]+)\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)\s*'
        '([0-9\.e\-\+]+)\s*([0-9\.e\-\+]+)',
        re.MULTILINE)
    with open(datafile,"r") as errors:
        errors = errors.read()
        for (ls, ks, la, ka, n, l2error, aposteriori) \
                in dataPattern.findall(errors):
            ls = int(ls)
            ks = int(ks)
            la = int(la)
            ka = int(ka)
            n  = int(n)
            l2error     = float(l2error)
            aposteriori = float(aposteriori)
            insertData(data, ls, ks, la, ka, n, l2error, aposteriori)

def plot(data,
         dataselect=lambda d: d['l2'],
         labelmaker=lambda ls, ks, la, ka: \
                           str(ls)+' '+str(ks)+' '+str(la)+' '+str(ka),
         outputfile='error.pdf',
         title=None,
         xlabel='mesh size $H$',
         ylabel='$L_2$ errors',
         xlim=None,
         ylim=None,
         xscale='log',
         yscale='log',
         legendlocation='upper left'):
    fig = plt.figure()
    if title != None:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    for ls in sorted(data):
        for ks in sorted(data[ls]):
            for la in sorted(data[ls][ks]):
                for ka in sorted(data[ls][ks][la]):
                    d = data[ls][ks][la][ka]
                    mesh = map(lambda x: 1./x, d['n'])

                    line = plt.plot(mesh, dataselect(d))
                    if isinstance(labelmaker, str) or labelmaker==None:
                        label=labelmaker
                    else:
                        label=labelmaker(ls,ks,la,ka)
                    plt.setp(line, linewidth=2.0,
                             marker='o', markersize=3.0,
                             label=label)

    plt.xscale(xscale)
    plt.yscale(yscale)
    if labelmaker != None:
        plt.legend(loc=legendlocation,shadow=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)

    plt.clf()


datadir = '.'
datafiles = map(lambda s: datadir + '/' + s,
                filter(lambda s: s.startswith('convergence_error_'),
                       os.listdir(datadir)))
data = dict()
for datafile in datafiles:
    readData(data, datafile)

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':16})
mpl.rc('text', usetex=True)

plot(data,
     dataselect=lambda d: d['l2'],
     labelmaker=lambda ls, ks, la, ka: \
                       '$\\mathbb{V}_h$, with $h=2^{-' + str(ls) + '}H$',
     outputfile='error.pdf',
     # title='$L_2$ error of $\\varphi$',
     ylabel='$L_2$ errors',
     ylim=[1e-4,1e-1])

plot(data,
     dataselect=lambda d: d['apost'],
     labelmaker=lambda ls, ks, la, ka: \
                       '$\\mathbb{V}_h$, with $h=2^{-' + str(ls) + '}H$',
     outputfile='error_aposteriori.pdf',
     # title='a posteriori error of $u=(\\varphi, \\theta)$',
     ylabel='a posteriori errors')

plot(data,
     dataselect=lambda d: map(lambda (l, a): a/l,
                              zip(d['l2'], d['apost'])),
     labelmaker=lambda ls, ks, la, ka: \
                       '$\\mathbb{V}_h$, with $h=2^{-' + str(ls) + '}H$',
     outputfile='error_apost_rel.pdf',
     # title='relative a posteriori error' \
     #       '$\\frac{e_{a posteriori}}{e_{exact}}$',
     ylabel='relative a posteriori errors',
     yscale='linear',
     legendlocation='lower left')
