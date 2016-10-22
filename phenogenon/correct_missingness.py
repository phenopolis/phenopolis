#!/bin/env python

'''
this is to remove biased missingness in variants in gene-hpo analysis

the method:
    1. walk through gene_hpo/release/version/EN*.json files
    2. read json, for each hpo, for each variant, check its location's missingness
    3. if missed on a given individual, check KINGSHIP imputed data, if available, impute. if not,  
    4. re-check if the patient is still comp-het on this gene if the data come from comp_het
    5. re-calculate fisher_p
'''
from optparse import OptionParser
from lookups import *
import pymongo
import json
import math
from scipy.stats import chi2_contingency
import numpy as np
import os
import sys
import time
import rest
import fisher
import copy
import pysam
from draw_hpo_graph import *

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db=conn['hpo']
db = conn['uclex']
patient_db=conn['patients']

def check_var(var,pat_h_array):
    # return pvalue
    missed_nohpo = called_nohpo = missed_hpo = called_hpo = 0
    chrom,pos,ref,alt = var.split('-')
    record = tbx.fetch(str(chrom),int(pos)-1,int(pos))
    for row in record:
        row=row.split('\t')
        # genotype starts from row[9]
        for ind in range(9,len(tbx_header)):
            if tbx_header[ind] not in pat_a_array: continue
            called = 0 if row[ind][0]==row[2]=='.' else 1
            in_hpo = 1 if tbx_header[ind] in pat_h_array else 0
            if (not called) and (not in_hpo):
                missed_nohpo += 1
            elif called and (not in_hpo):
                called_nohpo += 1
            elif (not called) and in_hpo:
                missed_hpo += 1
            elif called and in_hpo:
                called_hpo += 2
            else:
                raise 'something is wrong when trying to count miss call and hpo numbers'
    # do some stats
    obs = np.array([[missed_nohpo,called_nohpo], [missed_hpo, called_hpo]])
    # use chi square without Yates' correction to be conserved
    p_val = chi2_contingency(obs,correction=False)[1]
    return p_val


def cleanse(data,mode):
    # cleanse variants whose miss calls correlate with the hpo appearance
    # mark them as failed
    # note that BAD records are HPO specific!!
    data['style']='strict_'+style
    for hpo,value in data.iteritems():
        if hpo == 'HP:0000001': continue # no need to process ALL
        # get pat_h
        pat_h_array = hpo_freq[hpo]['patients']
        value['bad_vars'] = []
        checked = []
        # have to check all the patients in gene
        pat_g = pat_gh = 0
        pat_h = len(pat_h_array)
        for p,val in data['HP:0000001']['data'].iteritems():
            good_vars_count = 0
            for var in val['var']:
                if var['variant_id'] in value['bad_vars']:
                    var['failed'] = 1
                    continue
                if var['variant_id'] in checked:
                    continue
                # use chi square without Yates' correction to be conserved
                p_val = check_var(var['variant_id'],pat_h_array)
                if p_val <= 0.05:
                    # bad position, throw away
                    value['bad_vars'].append(var['variant_id'])
                else:
                    good_vars_count += 1
                checked.append(var['variant_id'])
            # is this patient good?
            if mode == 'hom_comp' and good_vars_count > 1:
                pat_g += 1
                if p in pat_h_array:
                    pat_gh += 1
            elif mode == 'het' and good_vars_count:
                pat_g += 1
                if p in pat_h_array:
                    pat_gh += 1
        # calculate adjusted p value 
        p_val = fisher.pvalue(pat_a-pat_h-pat_g+pat_gh,pat_h-pat_gh,pat_g-pat_gh,pat_gh)
        value['adj_p_val'] = p_val.right_tail
        value['adj_pat_g'] = pat_g
        value['adj_pat_gh'] = pat_gh
                
            
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)

parser.add_option("--chrom",
                  dest="chrom",
                  help="which chrom to process?")
(options, args) = parser.parse_args()
version = 2
stamp = '2016_Aug'
style = 'related'
vcf_location = '/SAN/vyplab/UCLex/current'
vcf_release = 'July2016'
#genes = get_chrom_genes([options.chrom], db)
genes = 'ABCA4'
hpo_freq = get_hpo_size_freq('hpo_freq_'+style+'_'+stamp+'_2.tsv')
pat_a_array = hpo_freq['HP:0000001']['patients']
pat_a = len(pat_a_array)
vcf_file = os.path.join(vcf_location,'mainset_'+vcf_release+'_chr'+options.chrom+'.vcf.gz')
tbx = pysam.TabixFile(vcf_file) 
# get tbx header
tbx_header = []
for h in tbx.header:
    if h[:2] == '##': continue
    tbx_header = h.split('\t')
    break

#genes=['HGSNAT'];
i = 0
for gene_id in genes:
    i += 1
    
    if not gene_id.startswith('ENSG'): gene_id = get_gene_by_name(db, gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    # already done?
    filename = os.path.join('gene_hpo',stamp,'strict_'+style,str(version),gene_id + '.json')
    if os.path.isfile(filename) and os.stat(filename).st_size:
        bad = 0
        try:
            inf = open(filename,'r')
            test_json = json.load(inf)
            inf.close()
        except ValueError:
            bad = 1
            # truncated json
        if not bad and test_json.get('version', None) == version:
            continue

    print '====='
    print gene_id
    print '%s / %s' % (i, len(genes))
    infile = os.path.join('gene_hpo',stamp,'loose_'+style,str(version),gene_id+'.json')
    data = json.load(open(infile,'r'))
    for mode in ['hom_comp','het']:
        if not data[mode]: continue
        cleanse(data[mode],mode)
        #filename = os.path.join(gene_id + '_' + mode + '.json')
    outf = open(filename,'w')

