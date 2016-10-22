#!/bin/env python
'''
this script is to snapshot patient hpo information for further analyses, such as gene-hpo relationship
'''
import sys
import pymongo
from lookups import *

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db=conn['hpo']
db = conn['uclex']
patient_db=conn['patients']

release = '2016_Aug'
outf = open('patients_hpo_snapshot_'+release+'.tsv','w')


unrelated = open('/SAN/vyplab/UCLex/mainset_July2016/kinship/UCL-exome_unrelated.txt','r').readlines()
unrelated = [i.rstrip() for i in unrelated]

outf.write('#p_id   unrelated   HPO\n')

for p in patient_db.patients.find({},{'features':1, 'external_id':1}):
    #p_id   unrelated   hpo
    if 'features' not in p:
        continue
    unrelated_flag = 1 if p['external_id'] in unrelated else 0
        
    hpos = [i['id'] for i in p['features'] if i['observed'] == 'yes']
    if not hpos:
        continue
    # replace obsolete hpos
    hpos = [replace_hpo(hpo_db, [h,h])[0] for h in hpos]
    all_hpos = [] # union of all ancestors
    for hpo in hpos:
        anc = get_hpo_ancestors(hpo_db, hpo)
        for a in anc:
            all_hpos.extend(a['id'])
    all_hpos = set(all_hpos)
    #write to file
    outf.write('%s\t%s\t%s\n' % (p['external_id'],unrelated_flag,','.join(all_hpos)))
