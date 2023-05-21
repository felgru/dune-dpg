#!/usr/bin/python3
import argparse

from output_parser import readData

aparser = argparse.ArgumentParser(
        description='print a Latex table of the given convergence results')
aparser.add_argument('--standalone',
                     action='store_true', default=False,
                     help='print standalone document that can be directly '
                          'compiled by pdflatex')
aparser.add_argument('infile', action='store')
args = aparser.parse_args()

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
          r'\newcommand{\cT}{\ensuremath{\mathcal{T}}}'
          '\n\n'
          r'\begin{document}')

def print_table(data):
    print(r'\begin{tabular}{r|l|r}')
    print(r'iteration & a posteriori error & \#DoFs \\')
    print(r'\hline')
    print_rows(data)
    print(r'\end{tabular}' '\n')

def print_rows(data):
    for n, err, dofs in zip(data['iterationIndices'],
                            data['globalAccAposteriori'],
                            data['dofs']):
        print(rf'{n} & {err} & {dofs} \\')

data = readData(args.infile)
if args.standalone:
    print_preamble()
print_table(data)
if args.standalone:
    print(r'\end{document}')
