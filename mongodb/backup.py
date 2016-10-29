#! /home/rmhanpo/miniconda2/bin/python
import sys

password_file=sys.argv[1]
from phenotips_python_client import PhenotipsClient
conn=PhenotipsClient(host='localhost',port=8080) 
password=file(password_file,'r').read().strip() 
conn.dump_to_mongodb(auth='Admin:%s'%password,mongo_host='localhost',mongo_port='27017',mongo_dbname='patients')


import sys
import pymongo

'''
Build an HPO to patient id cache, nightly, to facilitate rapid lookup of individuals by HPO term.
'''

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db=conn['hpo']
patients_db=conn['patients']

patients_db.hpo_cache.drop()
patients_db.hpo_cache.create_index('hpo_id',unique=True)

def get_hpo_patients(hpo_db, patients_db, hpo_id, cached=True,verbose=False):
    """
    Get patients with HPO term.
    """
    if cached:
        return [p for p in patients_db.patients.find({'external_id':{'$in':patients_db.hpo_cache.find_one({'hpo_id':hpo_id})['external_id']}}) if 'external_id' in p]
    if 'HP:0000001' == hpo_id: return [p for p in patients_db.patients.find() if 'external_id' in p]
    patients = [p for p in patients_db.patients.find({'features.id':hpo_id}) for f in p['features']  if f['id']== hpo_id and f['observed']=='yes']
    if verbose: print(hpo_id,len(patients))
    for r in hpo_db.hpo.find({'is_a':hpo_id}):
        for i in r['id']: patients+=list(itertools.chain(get_hpo_patients(hpo_db,patients_db,i,cached=cached,verbose=verbose))) 
    #remove duplicates
    patients={v['external_id']:v for v in patients}.values()
    return patients

for hpo in hpo_db.hpo.find():
    hpo_patients=get_hpo_patients(hpo_db, patients_db, hpo['id'][0], cached=False, verbose=False)
    print(patients_db.hpo_cache.insert( {'hpo_id':hpo['id'][0], 'external_id':[p['external_id'] for p in hpo_patients if 'external_id' in p]} ))



