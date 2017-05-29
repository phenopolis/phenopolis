from __future__ import print_function
import sys
import pymongo
import json
import re

host='phenotips.cs.ucl.ac.uk'

conn = pymongo.MongoClient(host=host, port=27017)
db=conn['patients']

headers=["Internal reference number or ID","Chromosome","Start","Genome assembly","Reference allele","Alternate allele","Transcript","Gene name","Intergenic","Chromosomal sex","Open-access consent","Age at last clinical assessment","Prenatal age in weeks","Note","Inheritance","Pathogenicity","Phenotypes","HGVS code","Genotype","Responsible contact"]
#headers=["Internal reference number or ID","Genome assembly","Gene name","Open-access consent","Phenotypes","Responsible contact"]

print(','.join(map(lambda x: '"%s"'%x,headers)))
for p in db.patients.find({'external_id':{'$regex':re.compile('IRDC_.*_LON_.*')}}):
    r=dict()
    if 'GC' not in p['external_id']: continue
    if 'genes' in p:
        r["Gene name"]= ', '.join([g['gene'] for g in p['genes']])
    else:
        r["Gene name"]= ''
    r["Internal reference number or ID"]=re.match('.*_(GC.*)',p['external_id']).group(1)
    r["Chromosome"]=''
    r["Start"]=''
    r["Genome assembly"]='GRCh37/hg19'
    r["Reference allele"]=''
    r["Alternate allele"]=''
    r["Transcript"]=''
    r["Intergenic"]=''
    r["Chromosomal sex"]=''
    r["Open-access consent"]='No'
    r["Age at last clinical assessment"]=''
    r["Prenatal age in weeks"]=''
    r["Note"]=''
    r["Inheritance"]=''
    r["Pathogenicity"]=''
    r["Phenotypes"]=', '.join([f['id'] for f in p['features'] if f['observed']=='yes'])
    r["HGVS code"]=''
    r["Genotype"]=''
    r["Responsible contact"]='Andrew Webster'
    print(','.join(['"%s"' % r[k] for k in headers]))

