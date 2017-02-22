from os import listdir, chdir
from os.path import isfile, join
import pymongo

conn = pymongo.MongoClient(host='localhost', port=27017)
db=conn['patients']

print ','.join(['eid', 'gender', 'observed_features', 'non_observed_features', 'genes', 'solved'])
for p in db.patients.find():
    if 'external_id' not in p: continue
    eid=p['external_id']
    observed_features=[f['id'] for f in p.get('features',[]) if f['observed']=='yes']
    non_observed_features=[f['id'] for f in p.get('features',[]) if f['observed']=='no']
    gender = p['sex']
    family_history = p.get('family_history',[])
    global_mode_of_inheritance = [i['id'] for i in p.get('global_mode_of_inheritance',[])]
    observed_features=list(set(observed_features+global_mode_of_inheritance))
    genes = [g['gene'] for g in p.get('genes',[])]
    solved = p.get('solved',dict()).get('status','')
    print ','.join([eid, gender, ';'.join(observed_features), ';'.join(non_observed_features), ';'.join(genes), solved])
