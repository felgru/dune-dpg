#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import argparse
import sys, os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def createArrowPlots(datafile):
    dirPattern = re.compile(
        r'Directions are:\n'
        r'((\-?[0-9]+\.?[0-9]*e?-?[0-9]* \-?[0-9]+\.?[0-9]*e?-?[0-9]*\n)*)'
        , re.MULTILINE)
    with open(datafile,"r") as errors:
        errors = errors.read()
        directions = dirPattern.findall(errors)
        # Number of iterations
        N = len(directions)
        # List of directions (for all iterations)
        S = []
        for n in range(N):
            # List of directions for iteration n
            sIter = []
            dirSplit = (directions[n][0]).split()
            for idir in range(int(len(dirSplit)/2)):
                s = { 'angle': [float(dirSplit[2*idir]), float(dirSplit[2*idir+1])] , 'index': idir}
                sIter.append(s)
                drawArrowPlot(n,idir,float(dirSplit[2*idir]),float(dirSplit[2*idir+1]))
            S.append(sIter)
    return [N, S]

def drawArrowPlot(n, dirNum, lx, ly):
    x0 = 0.
    y0 = 0.
    plt.quiver(x0, y0, lx, ly, angles='xy', scale_units='xy', scale = 1)
    plt.axis([-1, 1, -1, 1])
    plt.savefig('dir-n'+str(n)+'-dir'+str(dirNum)+'.pdf')
    plt.clf()

def reorder(S):
    reorderedS = []
    for n in range(len(S)):
        s_to_south = [s for s in S[n] if s['angle'][1]<=0.]
        s_to_north = [s for s in S[n] if s['angle'][1]>0.]
        reorder_s_to_south = sorted(s_to_south, key=lambda x: x['angle'][0])
        reorder_s_to_north = sorted(s_to_north, key=lambda x: -x['angle'][0])
        reorderedS.append(reorder_s_to_south+reorder_s_to_north)
    return reorderedS

def mergePlots(N, reorderedS
            , convert_u_png_to_pdf = False
            , prefixU = 'u-elevate'
            , prefixDir = 'dir'
            , prefixMerge = 'merge-elevate'
            , prefixSummary = 'summary-elevate'):
    # We order files depending on directions
    for n in range(N):
        for s in reorderedS[n]:
            # Old and new indexes for current direction
            oldindex = s['index']
            newindex = reorderedS[n].index(s)
            # Definition of file names
            zerosold = ''
            zerosnew = ''
            if oldindex<10:
                zerosold = '000'
            else:
                if oldindex<100:
                    zerosold = '00'
                else:
                    zerosold = '0'
            if newindex<10:
                zerosnew = '000'
            else:
                if newindex<100:
                    zerosnew = '00'
                else:
                    zerosnew = '0'
            # File names for u and direction (they have the old indexes)
            filenameU   = prefixU+'-n'+str(n)+'-dir'+zerosold+str(oldindex)
            filenameDir = prefixDir+'-n'+str(n)+'-dir'+str(oldindex)
            # File name for merged file (it has the new index)
            fileMerge   = prefixMerge+'-n'+str(n)+'-dir'+zerosnew+str(newindex)

            # We convert png file of u into a pdf
            # (because Visit cannot save in pdf)
            if convert_u_png_to_pdf:
                cmd = 'convert '+filenameU+'.png '+filenameU+'.pdf'
                os.system(cmd)
            # We merge the two files
            cmd = 'pdfjam '+filenameU+'.pdf'+' '+filenameDir+'.pdf -o '+fileMerge+'.pdf --nup 2x1 --twoside --landscape --scale 1.1 --delta "-2.0cm 0.cm" --offset \'2cm 0cm\''
            os.system(cmd)
        # We put together the merged files of each iteration
        cmd = 'pdfjam '+prefixMerge+'-n'+str(n)+'-dir*.pdf -o '+prefixSummary+'-n'+str(n)+'.pdf --landscape'
        os.system(cmd)

def makeVideo(N
    , prefixSummary = 'summary-elevate'
    , prefixVideo = 'video-elevate'):
    for n in range(N):
        # Make a gif movie
        cmd = 'convert -verbose -delay 100 -loop 1 -density 300 '+prefixSummary+'-n'+str(n)+'.pdf '+prefixVideo+'-n'+str(n)+'.gif'
        os.system(cmd)

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Parser
description = 'Creates plots and videos to visualize solution. \
Files u-n*-dir*.png or u-elevate-n*-dir*.png from Visit required.'
aparser = argparse.ArgumentParser(description=description)
# Arguments parser
aparser.add_argument("--elevate", type=str2bool, nargs='?',
                        const=True, default=True,
                        help="Work with elevated pictures.")
aparser.add_argument('infile', action='store', help='output file created by running periter code')
args = aparser.parse_args()

# Prefix depending on elevate mode
prefixU = ''
prefixDir = ''
prefixMerge = ''
prefixSummary = ''
prefixVideo = ''

if args.elevate:
    prefixU = 'u-elevate'
    prefixDir = 'dir'
    prefixMerge = 'merge-elevate'
    prefixSummary = 'summary-elevate'
    prefixVideo = 'video-elevate'
    # Info print
    print("Creating visualization images (elevate mode on)")
else:
    prefixU = 'u'
    prefixDir = 'dir'
    prefixMerge = 'merge'
    prefixSummary = 'summary'
    prefixVideo = 'video'
    # Info print
    print("Creating visualization images (elevate mode off)")

# Create plots with directions
[N, S] = createArrowPlots(args.infile)
# Merge plots with the direction and its solution u
reorderedS = reorder(S)
mergePlots(N
    , reorderedS
    , convert_u_png_to_pdf=False
    , prefixU=prefixU
    , prefixDir=prefixDir
    , prefixMerge=prefixMerge
    , prefixSummary=prefixSummary
    )
# Make a video
makeVideo(N
    , prefixSummary=prefixSummary
    , prefixVideo=prefixVideo)
