#!/bin/env python
'''
extract frequency from patients.dot file
will be replaced with a proper function.
output: 
HP:0000001  1.
HP:0000002  0.7442
'''
import os
import re

infile = os.path.join('static','dot','patients.dot')
outfile = 'hpo_freq.tsv'

inf = open(infile, 'r')
outf = open(outfile, 'w')

freq = None
for l in inf:
    if re.match('digraph', l) or re.search('->', l):
        continue
    
    m = re.match(r'(\d+) / (\d+)', l)
    if m:
        freq = float(m.groups()[0]) / float(m.groups()[1])
    else:
        if freq != None:
            hid = re.match(r'(HP:\d+)"', l).groups()[0]
            outf.write('%(hid)s\t%(freq)s\n' % locals())
            freq = None
