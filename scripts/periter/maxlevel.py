#!/usr/bin/python

from __future__ import print_function
import re
import sys

if len(sys.argv) != 2:
    print("You have to provide exactly one file as input!")
    sys.exit(1)

maxlevel = 0
with open(sys.argv[1],'r') as f:
    output = f.read()
    for m in re.finditer(r'  - Grid level: ([0-9]+)', output):
        level = int(m.group(1))
        maxlevel = max(level, maxlevel)
print("max level =", maxlevel)
