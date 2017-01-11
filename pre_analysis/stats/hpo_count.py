
from __future__ import print_function
import sys
from collections import Counter
import pymongo
import itertools

conn = pymongo.MongoClient(host='localhost', port=27017)
patients_db=conn['patients']
hpo_db=conn['hpo']


#
def get_hpo_patients(hpo_id):
    patients=[p for p in patients_db.patients.find({'features.id':hpo_id})]
    for r in hpo_db.hpo.find({'is_a':hpo_id}):
        for i in r['id']: patients+=list(itertools.chain(get_hpo_patients(i))) 
    #remove duplicates
    patients={v['external_id']:v for v in patients}.values()
    return patients

hpo_patients=dict()
for hpo_id in [hpo['id'] for hpo in hpo_db.hpo.find({'is_a':'HP:0000118'},{'id':1,'_id':0})]:
    hpo_id=hpo_id[0]
    hpo_patients[hpo_id]=len([p['external_id'] for p in get_hpo_patients(hpo_id)])

for hpo_id in hpo_patients:
    hpo_term=hpo_db.hpo.find_one({'id':hpo_id},{'name':1,'_id':0})['name'][0]
    print(hpo_term, hpo_id, hpo_patients[hpo_id], sep=',')

