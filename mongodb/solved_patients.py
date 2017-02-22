import sys
import pymongo
import json

'''
Build an HPO to patient id cache to facilitate rapid lookup of individuals by HPO term.
'''


conn = pymongo.MongoClient(host='phenotips', port=27017)
db = conn['uclex']
patient_db=conn['patients']

db.solved_patients.drop()

db.solved_patients.create_index('external_id',unique=True)
db.solved_patients.create_index('genes')

for p in patient_db.patients.find():
    if 'external_id' not in p: continue
    #print p['external_id'], p['solved'], p.get('genes',[])
    if 'genes' not in p: continue
    solved_genes=dict()
    for g in p['genes']:
        het=[]
        hom=[]
        for var in db.variants.find({'canonical_gene_name_upper':g['gene']}):
            if p['external_id'] in var['het_samples'] and var['HET_COUNT']<20:
                het+=[var['variant_id']]
                print(p['external_id'], g['gene'], 'het',var['variant_id'],var['most_severe_consequence'],var['HET_COUNT'],var['canonical_cadd'],var.get('kaviar',''))
            if p['external_id'] in var['hom_samples'] and var['HET_COUNT']<20:
                hom+=[var['variant_id']]
                print(p['external_id'], g['gene'],'hom',var['variant_id'],var['most_severe_consequence'],var['HET_COUNT'],var['canonical_cadd'],var.get('kaviar',''))
            solved_genes[g['gene']]={'het':het, 'hom':hom}
    print(db.solved_patients.insert({'external_id':p['external_id'], 'genes':solved_genes}))

