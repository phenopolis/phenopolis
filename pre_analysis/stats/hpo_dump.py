
from __future__ import print_function
import sys
from collections import Counter
import pymongo

conn = pymongo.MongoClient(host='localhost', port=27017)
db=conn['uclex-old']

print('external_id','id','label',sep='\t')
for p in db.patients.find():
    if 'features' in p:
        for f in p['features']:
            if f['observed']=='yes':
                print(p['external_id'],f['id'],f['label'],sep='\t')


