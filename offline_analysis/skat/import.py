from __future__ import print_function
import sys
import re
import argparse
import pymongo
import os
import csv
import glob


parser=argparse.ArgumentParser(description='Arguments to import.py')
parser.add_argument('--hpo', required=True)
parser.add_argument('--basedir', required=True)
args=parser.parse_args()

hpo=args.hpo
basedir=args.basedir

conn = pymongo.MongoClient()
db=conn['uclex']
#db.skat.remove()

def get_files(hpo,mode):
    hpo_dir=os.path.join(basedir,hpo)
    pattern=hpo_dir+'/'+hpo+'_chr_*_'+mode
    print(pattern)
    for d in glob.glob(pattern):
        if os.path.isfile(os.path.join(hpo_dir,d,'skat.csv')): yield os.path.join(hpo_dir,d,'skat.csv')
        else: continue

for mode in ['rec','dom']:
    for f in get_files(hpo,mode):
        d=csv.DictReader(file(f,'r'))
        for r in d:
            if r['FisherPvalue']=='NA': continue
            if float(r['SKATO'])>0.5: continue
            mongo_record={}
            mongo_record.update(r)
            mongo_record['variants']=[v.replace('_','-') for v in r['SNPs'].split(';')]
            del mongo_record['SNPs']
            mongo_record['HPO']=os.path.basename(hpo)
            for k in ['FisherPvalue','SKATO','MeanCallRateCtrls','MeanCallRateCtrls','minCadd','OddsRatio']:
                try:
                    if mongo_record[k]=='NA': mongo_record[k]='nan'
                    mongo_record[k]=float(mongo_record[k])
                except:
                    continue
            mongo_record=dict([(k.replace('.','_'),mongo_record[k]) for k in mongo_record])
            mongo_record['mode']=mode
            print(mongo_record)
            #for k in mongo_record: print(k)
            db.skat.insert(mongo_record)


db.skat.create_index('ENSEMBL')
db.skat.create_index('CompoundHetPvalue')
db.skat.create_index('End')
db.skat.create_index('HWEp')
db.skat.create_index('Start')
db.skat.create_index('Chr')
db.skat.create_index('min_depth')
db.skat.create_index('nb_alleles_cases')
db.skat.create_index('case_maf')
db.skat.create_index('Symbol')
db.skat.create_index('nb_ctrl_homs')
db.skat.create_index('nb_case_homs')
db.skat.create_index('MaxMissRate')
db.skat.create_index('nb_alleles_ctrls')
db.skat.create_index('nb_snps')
db.skat.create_index('nb_cases')
db.skat.create_index('minCadd')
db.skat.create_index('MeanCallRateCtrls')
db.skat.create_index('MeanCallRateCases')
db.skat.create_index('OddsRatio')
db.skat.create_index('MinSNPs')
db.skat.create_index('FisherPvalue')
db.skat.create_index('HPO')
db.skat.create_index('nb_ctrl_hets')
db.skat.create_index('variants')
db.skat.create_index('total_maf')
db.skat.create_index('MaxCtrlMAF')
db.skat.create_index('ctrl_maf')
db.skat.create_index('nb_ctrls')
db.skat.create_index('nb_case_hets')
db.skat.create_index('SKATO')
db.skat.create_index('maxExac')







