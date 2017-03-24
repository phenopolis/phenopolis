

from os import listdir, chdir
from os.path import isfile, join
import pymongo

class Individual(object):
    def __init__(self, filename, db=None, hpo='HP:0000001'):
        pass
    def load_individual(self):
        conn = pymongo.MongoClient(host='localhost', port=27017)
        db=conn['uclex-old']
        for p in db.patients.find():
            eid=p['external_id']
            print eid
            for k in ['all_variants', 'rare_variants', 'compound_hets', 'homozygous_variants']:
                if k in p:
                    variants_count=len(p[k])
                    print k, variants_count
                    print db.patients.update({'external_id':eid},{'$set':{k+'_count':variants_count}},upsert=True)
    
    def get_patient_observed_hpo(self):
        # returns [('HP:0000001', 'hell yeah')]
        this_patient = patient_db.patients.find_one({'external_id':patient}) 
        result = [(None, None)]
        if not this_patient:
            #print 'ERROR: %s not in patients db' % patient
            pass
        else:
            if 'features' not in this_patient:
                print 'WARNING: features not in ' + patient
            p_features = this_patient.get('features', [{'id':'HP:0000001', 'label':'All', 'observed': 'yes' }])
            result = [(f['id'], f['label']) for f in p_features if f['observed']=='yes']
        return result



