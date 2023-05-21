#!/usr/bin/python3
import argparse
import re
import subprocess

import matplotlib.pyplot as plt

def createArrowPlots(datafile):
    dirPattern = re.compile(
        r'Directions are:\n'
        r'((\-?[0-9]+\.?[0-9]*e?-?[0-9]* \-?[0-9]+\.?[0-9]*e?-?[0-9]*\n)*)'
        , re.MULTILINE)
    with open(datafile, "r") as errors:
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
    plt.savefig(f'dir-n{n}-dir{dirNum:0>4}.pdf')
    plt.clf()

def reorder(S):
    reorderedS = []
    for Sn in S:
        s_to_south = [s for s in Sn if s['angle'][1]<=0.]
        s_to_north = [s for s in Sn if s['angle'][1]>0.]
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
            # File names for u and direction (they have the old indexes)
            filenameU   = f'{prefixU}-n{n}-dir{oldindex:0>4}'
            filenameDir = f'{prefixDir}-n{n}-dir{oldindex:0>4}'
            # File name for merged file (it has the new index)
            fileMerge   = f'{prefixMerge}-n{n}-dir{newindex:0>4}'

            # We convert png file of u into a pdf
            # (because Visit cannot save in pdf)
            if convert_u_png_to_pdf:
                subprocess.run(['convert',
                                f'{filenameU}.png', f'{filenameU.pdf}'])
            # We merge the two files
            subprocess.run(['pdfjam', f'{filenameU}.pdf', f'{filenameDir}.pdf',
                            '-o', f'{fileMerge}.pdf',
                            '--nup', '2x1',
                            '--twoside', '--landscape',
                            '--scale', '1.1',
                            '--delta', '"-2.0cm 0.cm"',
                            '--offset', '"2cm 0cm"'])
        # We put together the merged files of each iteration
        subprocess.run(['pdfjam',
                        *sorted(Path('.').glob(f'{prefixMerge}-n{n}-dir*.pdf')),
                        '-o', f'{prefixSummary}-n{n}.pdf',
                        '--landscape'])

def makeVideo(N
    , prefixSummary = 'summary-elevate'
    , prefixVideo = 'video-elevate'):
    for n in range(N):
        # Make a gif movie
        subprocess.run(['convert', '-verbose',
                        '-delay', '100', '-loop', '1', '-density', '300',
                        f'{prefixSummary}-n{n}.pdf',
                        f'{prefixVideo}-n{n}.gif'])

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
aparser.add_argument('infile', action='store', help='output file created by running asti code')
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
