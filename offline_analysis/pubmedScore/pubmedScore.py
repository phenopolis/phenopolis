#!/bin/env python
'''
currently fetch_from_pubmed, parse_authors and parse_abstract are not in use
'''
from pprint import pprint
import os
import json
from Bio import Entrez
import pymongo
import sys
THIS_FILE_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(THIS_FILE_PATH,'..','commons'))
from phenopolis_utils import *
from time import gmtime, strftime
import time
import re
from optparse import OptionParser
import itertools
from urllib2 import HTTPError, URLError
import csv
import sys
import re
from collections import defaultdict, Counter

# get dbs
dbs = get_mongo_collections()
# username used to connect to pubmed
email = OFFLINE_CONFIG['pubmedScore']['email']

# Some constant
life = int(OFFLINE_CONFIG['pubmedScore']['lifetime']) #3600 * 24 * 30 for a month, 1 for fresh insert
now = time.mktime(time.localtime()) #get time in seconds

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
find the PID, Title or Abstract no matter what!
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

'''
query pubmed, and get result
'''
def pubmed_query(gene,keywords,lag=0,email='me@example.com'):
    Entrez.email = email
    # get reg
    reg = '\\b|\\b'.join(keywords)
    reg = '\\b' + reg + '\\b'
    reg = re.compile(reg, re.IGNORECASE)
    # get term
    term = gene + ' AND (' + ' OR '.join(['"' + t + '"' + '[Title/Abstract]' for t in keywords]) + ')'
    if lag:
        lag = int(lag/3600/24) # convert it to days, at least 1 day
        lag = 1 if lag < 1 else lag
            # need to update
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed', term=term, reldate=lag, datetype='pdat', usehistory='y'))
                break
            except URLError as err:
                print ('!!URLError %s' % err)
                time.sleep(2)
                attempt += 1
    else:
        # just search
        attempt = 1
        while attempt <= 10:
            try:
                search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=term, usehistory='y'))
                break
            except URLError as err:
                print ('URLError: %s at line 249' % err)
                time.sleep(2)
                attempt += 1
            except RuntimeError as err:
                print ('Runtime error: %s at line 276' % err)
                time.sleep(2)
                attempt += 1
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'score':0}
    # get search content
    if count:
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
            except SocketError as err:
                print('Socket error')
                time.sleep(2)
            except URLError as err:
                print ('URLError')
                time.sleep(2)
        record = Entrez.read(handle)
        if record:
            # got something. let's do some calculation
            for r in record['PubmedArticle']:
                # calculate score
                score = 0
                pid = str(find_item(r, 'PMID'))
                abstract_list = find_item(r, 'AbstractText')
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
                if title:
                    score = score + len(reg.findall(title))
                if abstract:
                    score = score + len(reg.findall(abstract))

                # add result to genes[gene_name]
                if score:
                    results['results'].append({
                        'id': pid,
                        'title': title,
                        'abstract': abstract,
                        'score': score
                    })
                    results['score'] = results['score'] + score
        results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results
'''
Extraction, pubmed search, and database insertion
dominant only used in the scoring

'''
def pubmed(gene_name, keywords, now, test=False):
    keywords.sort()
    # constructing searching term and cache database key
    term = '_'.join([gene_name.upper(), '-'.join(keywords).lower()])
    # check if the pubmed result already in the database, and if it is outdated
    saved = dbs['pubmedbatch'].cache.find_one({'key': term})

    lag = 0
    this_life = 1 if test else life
    if saved:
        lag = now - saved['date']
        # has record. let's see if it is out of date
        if lag  <= this_life:
            # up to date
            #print 'already in the database, and up to date'
            return saved
    #print 'number of results', len([r for r in results['results']])
    # update the database, and maybe results
    # update database now
    results = pubmed_query(gene_name,keywords,lag,email)
    if saved:
        results['results'].extend(saved['data'])
    if not test:
        dbs['pubmedbatch'].cache.update({'key': term}, {'$set': {'score':results['score'],'data':results['results'],'date':now}},upsert=True)
    return results


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
        pubmed(options.gene,dbs,keywords,now)
    elif options.patient:
        print options.patient
        patient=dbs['uclex'].patients.find_one({'external_id':options.patient})
        for v in patient['rare_variants']+patient['homozygous_variants']+patient['compound_hets']:
            pubmed(v['canonical_gene_name_upper'],dbs,keywords,now)
            keywords.sort()
            term = ','.join(keywords).lower()
            print(term)
            dbs['uclex'].patients.update({'external_id':options.patient},{'$set':{'pubmed_key':term}})
    else:
        patients = open(options.input, 'r')
        for p in patients:
            p = p.rstrip()
            print p
            for v in dbs['uclex'].patients.find_one({'external_id':p})['compound_hets']:
                pubmed(v['canonical_gene_name_upper'],dbs,keywords,now)
                keywords.sort()
                term = '-'.join(keywords).lower()
                dbs['uclex'].patients.update({'external_id':options.patient},{'$set':{'pubmed_key':term}})
