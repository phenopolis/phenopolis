#!/bin/env python
'''
connect to patients db, get hpo frequencies, write to hpo_freq.tsv
only consider unrelated individuals calculated by KING
'''
from lookups import *
import pymongo

unrelated = open('/slms/gee/research/vyplab/UCLex/mainset_July2016/kinship/UCL-exome_unrelated.txt','r').readlines()
unrelated = [i.rstrip() for i in unrelated]
outf = open('hpo_freqi.tsv', 'w')

def get_db(dbname=None):
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    return connect_db(dbname)

def connect_db(dbname=None):
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host='localhost', port=27017)
    print(client)
    if not dbname: dbname='uclex-old'
    print(dbname)
    return client[dbname]

hpo_db=get_db('hpo')
patient_db=get_db('patients')

result = {}
Sum = 0
for p in patient_db.patients.find({},{'features':1, 'external_id':1}):
    if 'features' not in p:
        continue
    #if p['external_id'] not in unrelated:
    #    continue
    hpos = [i['id'] for i in p['features'] if i['observed'] == 'yes']
    if not hpos:
        continue
    # replace obsolete hpos
    hpos = [replace_hpo(hpo_db, [h,h])[0] for h in hpos]
    Sum += 1
    all_hpos = [] # union of all ancestors
    for hpo in hpos:
        anc = get_hpo_ancestors(hpo_db, hpo)
        for a in anc:
            all_hpos.extend(a['id'])
    all_hpos = set(all_hpos)
    for hpo in all_hpos:
        if hpo in result:
            result[hpo]['patients'].append(p['external_id'])
        else:
            result[hpo]={'patients':[p['external_id']]}

print Sum
# get freq
for h, v in result.iteritems():
    l = len(set(v['patients']))
    #p = ','.join(v['patients'])
    outf.write('%(h)s\t%(l)s/%(Sum)s\n' % locals())
