#!/usr/bin/env python


'''
  Filter variants into patient collection directly from database.

  Runs exomiser prioritiser to score genes.

  Get genes with damaging variants (usually hom / het / rare) from individual
  and query them against pubmed using pubmedPython
'''

from pprint import pprint
import os
import json
import pymongo
import sys
import re
import itertools
from urllib2 import HTTPError, URLError
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import csv
from collections import defaultdict, Counter
#import rest as annotation
import vcf
from phenotips_python_client import PhenotipsClient
from optparse import OptionParser
import mygene
import lookups
import orm
import requests


def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item

def load_patient(individual,auth,hpo='HP:0000001',AC=6,kaviar=.05,consequence_exclude=['intron_variant','non_coding_transcript','5_prime_UTR_variant','3_prime_UTR_variant','upstream_gene_variant','downstream_gene_variant']):
    client = pymongo.MongoClient()
    hpo_db = client['hpo']
    db = client['uclex']
    patient_db = client['patients']
    patient_id=individual
    # Add patient to phenotips if it does not already exist
    pheno=PhenotipsClient()
    patient={u'features': {u'observed': u'yes', u'type': u'phenotype', u'id': hpo}, 'clinicalStatus': {u'clinicalStatus': u'affected'}, u'ethnicity': {u'maternal_ethnicity': [], u'paternal_ethnicity': []}, u'family_history': {}, u'disorders': [], u'life_status': u'alive', u'reporter': u'', u'genes': [], u'prenatal_perinatal_phenotype': {u'prenatal_phenotype': [], u'negative_prenatal_phenotype': []}, u'prenatal_perinatal_history': {u'twinNumber': u''}, u'sex': u'U', u'solved': {u'status': u'unsolved'}}
    eid=patient_id
    p=pheno.get_patient(auth=auth,eid=eid)
    if p: patient.update(p)
    else: raise 'patient not on phenotips'
    #patient_hpo_terms=','.join([f['id'] for f in patient['features'] if f['observed']=='yes'])
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, patient_id, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    patient['total_variant_count']=db.variants.find({'hom_samples':individual}).count()+db.variants.find({'het_samples':individual}).count()
    print('total_variant_count:',patient['total_variant_count'])
    def conditions(x):
        if 'AC' not in x or x['AC'] > AC: return False
        if x['EXAC'] and x['EXAC']['AC_POPMAX']>=AC: return False
        if 'kaviar' in x and x['kaviar']>kaviar: return False
        if 'canonical_gene_name_upper' not in x: return False
        if x['most_severe_consequence'] in consequence_exclude: return False
        return True
    def process(x):
        if type(x['canonical_gene_name_upper']) is list: x['canonical_gene_name_upper']=x['canonical_gene_name_upper'][0]
        gene_hpo_terms=lookups.get_gene_hpo(hpo_db,x['canonical_gene_name_upper'],False)
        gene_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in gene_hpo_terms])
        gene_hpo_ids=gene_hpo_terms.keys()
        #lookups.get_gene_hpo(hpo_db,gene_name,dot=False)
        #print 'gene', gene_hpo_ids
        #print 'patient', patient_hpo_ids
        common_hpo_ids=list(set(gene_hpo_ids) & set(patient_hpo_ids))
        # simplify hpo terms
        common_hpo_ids=lookups.hpo_minimum_set(hpo_db, common_hpo_ids)
        common_hpo_ids=[{'hpo_id':k,'hpo_term':patient_hpo_terms[k]['name']} for k in common_hpo_ids]
        x['HPO']=common_hpo_ids
        print x['canonical_gene_name_upper'],common_hpo_ids
        for g in x['genes']:
            r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
            if not r or not r['EntrezGene ID']: continue
            x['entrezgeneid']=r['EntrezGene ID']
            break
        if 'entrezgeneid' not in x:
            x['exomiser']={'humanScore':-1}
            return x
        url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        x['exomiser']=r.json()[0]
        print(x)
        return x
    patient['homozygous_variants']=[ process(x) for x in db.variants.find({'hom_samples':individual}) if conditions(x) ]
    patient['rare_homozygous_variants_count']=len(patient['homozygous_variants'])
    print('number of homozygous variants:',patient['rare_homozygous_variants_count'])
    patient['rare_variants']=[ process(x) for x in db.variants.find({'het_samples':individual}) if conditions(x) ]
    patient['rare_variants_count']=len(patient['rare_variants'])
    print('number of rare variants:',patient['rare_variants_count'])
    gene_counter=Counter([var['canonical_gene_name_upper'] for var in patient['rare_variants']])
    for var in patient['rare_variants']: var['gene_count']=gene_counter[var['canonical_gene_name_upper']]
    patient['compound_hets']=[process(x) for x in patient['rare_variants'] if x['gene_count']>1]
    patient['rare_compound_hets_count']=len(patient['compound_hets'])
    print('number of compound hets:',patient['rare_compound_hets_count'])
    patient["pubmedbatch_status"]=0
    patient["pubmed_key"]="blindness-macula-macular-pigmentosa-retina-retinal-retinitis-stargardt"
    db.patients.update({'external_id':patient_id}, patient, upsert=True)



if '__main__'==__name__:
    ###############################
    # parse options
    ###############################
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("--individual", dest="individual", help="id of patient in database")
    parser.add_option("--auth", dest="auth", help="login for Phenotips")
    parser.add_option("--hpo", default='HP:0000001', dest="hpo", help="HPO term if none found in Phenotips")
    (options, args) = parser.parse_args()
    individual=options.individual
    if not individual: parser.error('Individual needs to be specified')
    auth=options.auth
    hpo=options.hpo
    load_patient(individual,auth,hpo)


