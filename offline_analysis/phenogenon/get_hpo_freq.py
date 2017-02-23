#!/bin/env python
'''
connect to patients db, get hpo frequencies, write to hpo_freq.tsv
only consider unrelated individuals calculated by KING
'''
import json

release = '2017_Feb'
outfile = '../output/phenogenon/hpo_freq_%s.json' % release

infile = '../output/phenogenon/patients_hpo_snapshot_%s.tsv' % release
inf = open(infile, 'r')

result = {'related':{},'unrelated':{},'release':release}
for row in inf:
    if row[0] == '#': continue
    row = row.rstrip().split('\t')
    p = row[0]
    unrelated = int(row[1])
    hpos = row[2].split(',')
    for h in hpos:
        # union
        result['related'][h] = result['related'].get(h,[])
        result['related'][h].append(p)
        if unrelated:
            result['unrelated'][h] = result['unrelated'].get(h,[])
            result['unrelated'][h].append(p)

print 'related patients length: %s' % len(result['related']['HP:0000001'])
print 'unrelated patients length: %s' % len(result['unrelated']['HP:0000001'])

outf = open(outfile,'w')

json.dump(result,outf)
