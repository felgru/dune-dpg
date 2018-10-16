# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function)
from collections import defaultdict
import re

def readData(datafile):
    parametersPattern = re.compile(
        r'(PERITER|ASTI) algorithm\n'
        r'=================\n'
        r'Prescribed final accuracy: ([0-9]*\.?[0-9]*)\n'
        r'Henyey Greenstein kernel with gamma = ([0-9]*\.?[0-9]*)\n'
        r'Wavelet order: ([0-9]*\.?[0-9]*)\n'
        r'Kernel approximation with: (.*)\n'
        r'Maximum wavelet level: ([0-9]*\.?[0-9]*)\n'
        r'Maximum number of directions: ([0-9]*\.?[0-9]*)\n'
        r'(Periter|ASTI) parameters:\n'
        r'rho = ([0-9]*\.?[0-9]*)\n'
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
    etaPattern = re.compile(r'eta_n = \(1\+n\)\^{-alpha} rho\^n: ([0-9]*\.?[0-9]*e?[+-]?[0-9]+?)\n')
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
        parameters = { 'eps': parametersMatch.group(2)
                     , 'gamma':   parametersMatch.group(3)
                     , 'wltOrder':    parametersMatch.group(4)
                     , 'kernelApproxType':    parametersMatch.group(5)
                     , 'maxWltLevel':     parametersMatch.group(6)
                     , 'maxNumS': parametersMatch.group(7)
                     , 'rho': parametersMatch.group(9)
                     , 'kappa1': parametersMatch.group(10)
                     , 'kappa2': parametersMatch.group(11)
                     , 'kappa3': parametersMatch.group(12)
                     , 'CT': parametersMatch.group(13)
                     }
        # We always remove the first iteration step, as it is used
        # to initialize the scattering error estimate and thus has
        # scattering error = 0.
        singularValues = []
        if(parameters['kernelApproxType']=='SVD'):
            svPat = re.compile(r'([0-9]+\.?[0-9]*e?-?[0-9]*)\n', re.MULTILINE)
            singularValues = svPat.findall(singularValuesPattern.search(errors).group())[1:]
        iterationIndices = map(int, iterationIndicesPattern.findall(errors))[1:]
        etas = map(float, etaPattern.findall(errors))[1:]
        wltLevel = map(int, wltLevelPattern.findall(errors))[1:]
        numS = map(int, numSPattern.findall(errors))[1:]
        svdRank = map(int, svdRankPattern.findall(errors))[1:]
        matrixTH = matrixTHpattern.findall(errors)[1:]
        timeEvalKernel = map(float, timeEvalKernelPattern.findall(errors))[1:]
        aPost = map(float, aPostPattern.findall(errors))[1:]
        accKernel = map(float, accKernelPattern.findall(errors))[1:]
        globalAccIteratesDiff = map(float,
                globalAccIteratesDiffPattern.findall(errors))[1:]
        globalAccApriori = map(float, globalAccAprioriPattern.findall(errors))[1:]
        globalAccAposteriori = map(float,
                globalAccAposterioriPattern.findall(errors))[1:]
        dofs = map(int, dofsPattern.findall(errors))[1:]
        innerIterationsStats = defaultdict(list)
        for m in innerIterationsPattern.finditer(errors):
            innerIterationsStats[int(m.group(1))].append(
                    { 'numIterations': int(m.group(2))
                    , 'maxLevel': int(m.group(3))
                    , 'numDOFs': int(m.group(4))
                    })
        del innerIterationsStats[0]

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
           , 'globalAccIteratesDiff' : globalAccIteratesDiff
           , 'globalAccApriori': globalAccApriori
           , 'globalAccAposteriori': globalAccAposteriori
           , 'dofs': dofs
           , 'innerIterationsStats': innerIterationsStats
           }
