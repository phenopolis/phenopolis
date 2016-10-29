#!/bin/env python

from pprint import pprint
import os
import json
from Bio import Entrez
import pymongo
import sys
from time import gmtime, strftime
import time
import re
from optparse import OptionParser
import itertools
from urllib2 import HTTPError, URLError
import pysam
import csv
import sys
import re
from collections import defaultdict, Counter

Entrez.email = 'jing.yu@ndcn.ox.ac.uk'

# Some constant
today = strftime("%Y-%m-%d", gmtime())
life = 3600 * 24 * 30 #3600 * 24 * 30 for a month, 1 for fresh insert
now = time.mktime(time.localtime()) #get time in seconds
#life = now - 1470220103 + 600
HEADER = ['HUGO', 'HPO', 'consequence', 'ref(pubmedID)', 'description', 'OMIM', 'allele_freq', 'ExAC_freq', 'variant_id', 'p_change']
# get db
client = pymongo.MongoClient()
db_pubmed = client['pubmedbatch']
hpo_db = client['hpo']
db = client['uclex']
patient_db = client['patients']

###############################
# parse options
###############################

#input='patients.txt'
#patients = [options.patient] if options.patient else open(options.input, 'r')

#known_genes='ret_known_genes.txt',
#known_genes = open('ret_known_genes.txt', 'r').readline().strip().split()
#mask_genes = open('ret_blacklist.txt', 'r').readline().strip().split()


'''
update command line on what's going on
it actually does the remove
'''
def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()
    sys.stdout.write(' ' * 40)
    sys.stdout.flush()
    sys.stdout.write('\r')
    sys.stdout.flush()

'''
to check if an iterable is empty
'''
def peek(iterable):
    try:
        first = next(iterable)
    except RuntimeError:
        return None
    except StopIteration:
        return None
    return itertools.chain([first], iterable)


'''
find the freaking PID, Title or Abstract no matter what!
'''
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

"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def fetch_from_pubmed(query, lag=None):
    smashed_all=query['smashed_all']
    reg=query['reg']
    if lag:
        lag = lag/3600/24 # convert it to days
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed', term=smashed_all, reldate=lag, datetype='pdat', usehistory='y'))
                break
            except URLError as err:
                print ('PubmedBatch URLError: %s' % err)
                time.sleep(2)
                attempt += 1
    else:
        # just search
        attempt = 1
        while attempt <=10:
            try:
                print 'QUERY:', query
                search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=smashed_all, usehistory='y'))
                break
            except URLError as err:
                print ('URLError')
                time.sleep(2)
                attempt += 1
    # now done the search. let's get results
    count = int(search_results["Count"])
    results = {'results':[], 'total_score':0}
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed", restart=0, retmax=50, retmode="xml", webenv=search_results['WebEnv'], query_key=search_results['QueryKey'])
            break
        except HTTPError as err:
            if 500 <= err.code <= 599: print('Received error from server %s' % err)
            else: print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    if not peek(record): return []
    # got something. let's do some calculation
    for r in record:
        # calculate score
        score = 0
        pid = str(r['MedlineCitation']['PMID'])
        abstract_list = r['MedlineCitation']['Article']['Abstract']['AbstractText']
        # parse abstract
        abstract = ''
        if abstract_list:
            for a in abstract_list:
                if hasattr(a, 'attributes') and 'Label' in a.attributes:
                    abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                    abstract = abstract + a + '<br/>'
                else:
                    abstract = abstract + a
        title = find_item(r, 'ArticleTitle')
        if title: score = score + len(reg.findall(title))
        if abstract: score = score + len(reg.findall(abstract))
        # add result to genes[gene_name]
        if score:
            results['results'].append({
                'id': pid,
                'title': title,
                'abstract': abstract,
                'score': score
            })
            results['total_score'] = results['total_score'] + score
    results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    results['query']=query
    return results

'''
construct hpo_db to get HPO relative to eye
'''
def construct_hpo(individual):
    patient=db.patients.find_one({'external_id':individual})
    #print('number of variants',len(patient['variants']))
    #patient['variants']=db.variants.find({'variant_id':{'$in': map(lambda x: x.replace('_','-'), patient['variants']) }})
    #patient['variants']=map(lambda x: x.replace('_','-'), patient['variants'])
    patient['name']=patient['external_id']
    hpo_terms=[(f['id'],f['label'],) for f in patient_db.patients.find_one({'external_id':individual})['features'] if f['observed']=='yes']
    gene_hpo=dict()
    for hpo_id,hpo_term, in hpo_terms:
        #for gene_name in [x['gene_name'] for x in hpo_db.gene_hpo.find({'hpo_terms': hpo_id},{'gene_name':1,'_id':0})]:
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
    return gene_hpo


def parse_abstract(x):
    x['MedlineCitation']['Article']['Abstract']['AbstractText']


def parse_authors(x):
    for author in x['MedlineCitation']['Article']['AuthorList']: author

'''
Extraction, pubmed search, and database insertion
dominant only used in the scoring
'''
def pubmed(gene_name, keywords, now):
    keywords.sort()
    # make reg term for later use (calculating pubmed score)
    reg = '\\b|\\b'.join(keywords)
    reg = '\\b' + reg + '\\b'
    reg = re.compile(reg, re.IGNORECASE)
    # constructing searching term and cache database key
    smashed_OR = ['"' + t + '"' + '[Title/Abstract]' for t in keywords]
    smashed_OR = ' OR '.join(smashed_OR)
    smashed_term = ' AND (' + smashed_OR + ')'
    smashed_all = gene_name + smashed_term
    term = '_'.join([gene_name.upper(), '-'.join(keywords).lower()])
    # check if the pubmed result already in the database, and if it is outdated
    lag = 0 # use it as a flag of how to search. 0 = search; now-saved['date'] = update; 
    saved = db_pubmed.cache.find_one({'key': term})
    print(term)
    print(saved)
    if saved:
        lag = now - saved['date']
        # has record. let's see if it is out of date
        if lag  <= life:
            # up to date
            #print 'already in the database, and up to date'
            return None
    if lag > 0:
        lag = int(lag/3600/24) # convert to days
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed', term=smashed_all, reldate=lag, datetype='pdat', usehistory='y'))
                break
            except URLError as err:
                print ('PubmedBatch URLError: %s' % err)
                time.sleep(2)
                attempt += 1
    else:
        # just search
        attempt = 1
        while attempt <=10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=smashed_all, usehistory='y'))
                break
            except URLError as err:
                print ('URLError')
                time.sleep(2)
                attempt += 1
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed",
                                   restart=0,
                                   retmax=50,
                                   retmode="xml",
                                   webenv=search_results['WebEnv'],
                                   query_key=search_results['QueryKey']
                                   )
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print('Received error from server %s' % err)
            else:
                print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    count = int(search_results["Count"])
    #print 'count', count
    results = {'results':[], 'total_score':0}
    record = peek(record)
    print record
    if not record: 
        if saved:
            db_pubmed.cache.update({'key': term}, {'$set': {'date': now}})
        else:
            db_pubmed.cache.update({'key': term}, {'$set': {'score':results['total_score'],'data':results['results'],'date':now}},upsert=True)
        return None
    for r in record:
        # calculate score
        score = 0
        try:
            pid = str(r['MedlineCitation']['PMID'])
            abstract_list = r['MedlineCitation']['Article']['Abstract']['AbstractText']
        except:
            print 'problem with abstract list'
            continue
        # parse abstract
        abstract = ''
        if abstract_list:
            for a in abstract_list:
                if hasattr(a, 'attributes') and 'Label' in a.attributes:
                    abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                    abstract = abstract + a + '<br/>'
                else:
                    abstract = abstract + a
        title = find_item(r, 'ArticleTitle')
        if title: score = score + len(reg.findall(title))
        if abstract: score = score + len(reg.findall(abstract))
        # add result to genes[gene_name]
        if score:
            results['results'].append({
                'id': pid,
                'title': title,
                'abstract': abstract,
                'score': score,
                'vote':0, # customer opinion on the paper's relevance to the gene and the phenotypes
            })
            results['total_score'] = results['total_score'] + score
    #print 'number of results', len([r for r in results['results']])
    # update the database, and maybe results
    # update database now
    if saved:
        results['results'].extend(saved['data'])
    db_pubmed.cache.update({'key': term}, {'$set': {'score':results['total_score'],'data':results['results'],'date':now}},upsert=True)


if __name__ == '__main__':
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("-i", "--input",
                      dest="input", default='patients.txt',
                      help="patient list, new line separated: [default: %default]")

    parser.add_option("-g", "--gene",
                      dest="gene", default=None,
                      help="gene [default: %default]")

    parser.add_option('-p', '--patient',
                      dest='patient',
                      help='the patient id to update')

    parser.add_option('-k', '--keywords',
                      dest='keywords',
                      default="retina,retinal,retinitis,blindness,macula,macular,stargardt,pigmentosa",
                      help='')

    (options, args) = parser.parse_args()
    keywords=[i.strip() for i in options.keywords.split(',')]
    if options.gene:
        pubmed(options.gene,keywords,now)
    elif options.patient:
        print options.patient
        patient=db.patients.find_one({'external_id':options.patient})
        for v in patient['rare_variants']+patient['homozygous_variants']+patient['compound_hets']:
            pubmed(v['canonical_gene_name_upper'],keywords,now)
            keywords.sort()
            term = ','.join(keywords).lower()
            print(term)
            db.patients.update({'external_id':options.patient},{'$set':{'pubmed_key':term}})
    else:
        patients = open(options.input, 'r')
        for p in patients:
            p = p.rstrip()
            print p
            for v in db.patients.find_one({'external_id':p})['compound_hets']:
                pubmed(v['canonical_gene_name_upper'],keywords,now)
                keywords.sort()
                term = '-'.join(keywords).lower()
                db.patients.update({'external_id':options.patient},{'$set':{'pubmed_key':term}})
