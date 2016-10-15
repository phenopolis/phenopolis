#!/bin/env python
# calculate HPO-gene correlations, and rank
# p_hg = 1 - binom.cdf(pat_gh - 1, pat_g, r_ga)
# p_hg: the possibility that the correlation is random
# pat_gh: number of patients having both the genotype and the hpo terms
# pat_g: number of patients having the genotype
# r_ga: ratio of pat_g : pat_a (number of all the patients)
# all the parameters can be extracted by the HPO-gene json files

from optparse import OptionParser
import json
import math
from scipy.stats import binom
import os
import sys
from operator import itemgetter

for file in os.listdir("dot"):
    # for now i'm not filtering on cadd, since it is not up-to-date yet. Maybe will introduce it as a filter later in my life. do I have a long life? let's hope so.
    if file.endswith(".json") and file[:4] == 'ENSG':
        print(file)
        # getting gene id and mode
        match = file.split('_')
        gene = match[0]
        mode = 'hom_comp' if match[1] == 'hom' else 'het'
        if os.stat("dot/"+file).st_size == 0:
            continue
        # read in data
        inf = open('dot/'+file, 'r')
        data = json.load(inf)
        inf.close()
        r_ga = 0
        if 'data' in data:
            data = data['data']
        if not data:
            continue
        if 'data' not in data[data.keys()[0]]:
            continue
        pat_g = len(set([i for k,v in data.iteritems() for i in v['data'].keys()]))
        for hpo, value in data.iteritems():
            filename = hpo.split(':')[1]
            hpo_file = os.path.join('dot','hpo_gene',filename+'.json')
            # calculate r_ga
            if not r_ga:
                pat_a = float(value['expect_freq'].split('/')[1])
                r_ga = pat_g / pat_a

            # get pat_gh
            pat_gh = len(value['data'])
            pat_h = int(value['expect_freq'].split('/')[0])
            # calculate p value
            p_hg = 1 - binom.cdf(pat_gh - 1, pat_h, r_ga)
            # if p value > 0.05, forget about it
            if p_hg > 0.05:
                continue
            hpo_outf = open(hpo_file,'a')
            hpo_outf.write('%s\t%s\t%s\t%s/%s\t%s/%s\n' % (mode,gene,p_hg,pat_gh,pat_h,pat_g,pat_a))
            hpo_outf.close()

sys.exit()
for h,v in hpos.iteritems():
    # sort and write
    for mode in ['hom_comp','het']:
        # sort
        v['data'][mode] = sorted(v['data'][mode], key=itemgetter('p_value'))
    # file name doesn't like :, only use the number as file name
    filename = h.split(':')[1]
    outfile = os.path.join('dot','hpo_gene',filename+'.json')
    outf = open(outfile,'w')
    json.dump(v,outf)
