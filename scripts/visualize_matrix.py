#!/usr/bin/python

# This script can be used to generate a heatmap of a matrix
#
# Usage: visualize_matrix.py matrix_file
# where matrix_file is a text file that contains the output of the
# printmatrix function, e.g.
# printmatrix(of, stiffnessMatrix, "stiffnessMatrix", "--");
#
# The generated heatmap might be useful when debugging matrices to
# spot anomalies in the sparsity pattern.

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

input_name = sys.argv[1]
output_name = input_name + '.png'
with open(input_name,'r') as f:
    lines = f.readlines()
    dimensions = re.search(r'\[n=([0-9]+),m=([0-9]+),', lines[0])
    rows = int(dimensions.group(1))
    cols = int(dimensions.group(2))
    m = np.empty((rows, cols))
    for line in lines[1:]:
        line = line.split()
        row = int(line[1])
        col = 0
        for entry in line[2:]:
            if entry != '.':
                m[row,col] = np.float64(entry)
            else:
                m[row,col] = np.float64(0.)
            col += 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(m, norm=mpl.colors.LogNorm())
    fig.colorbar(cax)
    plt.savefig(output_name)
