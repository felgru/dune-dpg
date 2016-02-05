#!/usr/bin/python
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab

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
         title='$L_2$ error of $u$',
         xlabel='mesh size $h$',
         ylabel='$L_2$ errors',
         xscale='log',
         yscale='log',
         legendlocation='upper left'):
    fig = plt.figure()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

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
    plt.show()
    plt.savefig(outputfile)

    plt.clf()


datadir = '.'
datafiles = map(lambda s: datadir + '/' + s,
                filter(lambda s: s.startswith('convergence_error_'),
                       os.listdir(datadir)))
data = dict()
for datafile in datafiles:
    readData(data, datafile)

plot(data,
     dataselect=lambda d: d['l2'],
     labelmaker=lambda ls, ks, la, ka: \
                       'unrefined test space' if ls==0 \
                       else 'once refined test space' if ls==1 \
                       else 'twice refined test space' if ls==2 \
                       else str(ls)+' times refined test space',
     outputfile='error.pdf',
     title='$L_2$ error of $u$',
     ylabel='$L_2$ errors')

plot(data,
     dataselect=lambda d: d['apost'],
     labelmaker=lambda ls, ks, la, ka: \
                       'unrefined test space' if ls==0 \
                       else 'once refined test space' if ls==1 \
                       else 'twice refined test space' if ls==2 \
                       else str(ls)+' times refined test space',
     outputfile='error_aposteriori.pdf',
     title='aposteriori error of $u$, $\\theta$',
     ylabel='a posteriori errors')

plot(data,
     dataselect=lambda d: map(lambda (l, a): a/l,
                              zip(d['l2'], d['apost'])),
     labelmaker=lambda ls, ks, la, ka: \
                       'unrefined test space' if ls==0 \
                       else 'once refined test space' if ls==1 \
                       else 'twice refined test space' if ls==2 \
                       else str(ls)+' times refined test space',
     outputfile='error_apost_rel.pdf',
     title='relative a posteriori error $\\frac{e_{a posteriori}}{e_{exact}}$',
     ylabel='relative a posteriori errors',
     yscale='linear',
     legendlocation='lower left')
