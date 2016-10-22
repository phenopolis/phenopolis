#!/bin/env python
'''
re-calculate p using CADD cutoff 15
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
from plinkio import plinkfile

def write_to_data(data,dp,mode):
    # calculate and populate the p values
    # calculate pat g for related and unrelated cases
    related_pat = set([p for h,v in dp.iteritems() for p in v['data'] ])
    related_pat_g = len(related_pat)
    unrelated_pat_g = len([p for p in related_pat if p in patients_hpo['unrelated']]) 
    for h,v in dp.iteritems():
        related_pat_a = v['related_pat_a']
        unrelated_pat_a = v['unrelated_pat_a']
        related_pat_h = v['related_pat_h']
        unrelated_pat_h = v['unrelated_pat_h']
        related_pat_gh_set = set([p for p in v['data']])
        related_pat_gh = len(related_pat_gh_set)
        unrelated_pat_gh = len([p for p in related_pat_gh_set if p in patients_hpo['unrelated']])
        related_p_val = fisher.pvalue(related_pat_a-related_pat_h-related_pat_g+related_pat_gh,related_pat_h-related_pat_gh,related_pat_g-related_pat_gh,related_pat_gh)
        unrelated_p_val = fisher.pvalue(unrelated_pat_a-unrelated_pat_h-unrelated_pat_g+unrelated_pat_gh,unrelated_pat_h-unrelated_pat_gh,unrelated_pat_g-unrelated_pat_gh,unrelated_pat_gh)
        data[h]['related_'+mode+'_p_val'] = related_p_val.right_tail
        data[h]['unrelated_'+mode+'_p_val'] = unrelated_p_val.right_tail
        data[h]['related_'+mode+'_pat_g'] = related_pat_g
        data[h]['unrelated_'+mode+'_pat_g'] = unrelated_pat_g
        data[h]['related_'+mode+'_pat_gh'] = related_pat_gh
        data[h]['unrelated_'+mode+'_pat_gh'] = unrelated_pat_gh
def populate_fisher_p(data,mode):
    # add fisher p value into the data.
    # if mode == het, has to calculate d_p_val for dominant p
    exac_cutoff = 0.01
    if mode == 'het': exac_cutoff = 0.001

    if not data: return None
    dp = copy.deepcopy(data)
    for h,v in dp.iteritems():
        failed = []
        for p,var in v['data'].iteritems():
            count = 0
            for variant in var['var']:
                if variant['exac'] <= exac_cutoff and variant['cadd'] >= CADD_cutoff: count += 1
            if not count: failed.append(p)
        for p in failed:
            del v['data'][p]
    if mode == 'hom_comp':
        write_to_data(data,dp,'recessive')
    elif mode == 'het':
        # dominant_all doesn't eliminate hom_comphet
        write_to_data(data,dp,'dominant_all')
        # dominant_single, remove patients who have more than one var on gene
        for h,v in dp.iteritems():
            failed = []
            for p,var in v['data'].iteritems():
                count = 0
                for variant in var['var']:
                    if variant['exac'] <= exac_cutoff and variant['cadd'] >= CADD_cutoff: count += 1
                if count > 1:
                    failed.append(p)
            for p in failed:
                del v['data'][p]
        write_to_data(data,dp,'dominant_single')

def get_chrom_genes(chroms, db):
    # give chrom numbers, get all genes on them
    result = []
    for chrom in chroms:
        genes = [g['gene_id'] for g in db.genes.find({'chrom':str(chrom)})]
        result.extend(genes)
    return result

def parse_patients_hpo(file):
    inf = open(file,'r')
    result = {'related':{},'unrelated':{}}
    for row in inf:
        if row[0] == '#': continue
        row = row.rstrip().split('\t')
        result['related'][row[0]] = {'hpo':row[2].split(',')}
        if row[1] == '1':
            result['unrelated'][row[0]] = result['related'][row[0]]
    return result

def get_gradient_colour_raw(this_IC, max_IC, min_IC, max_colour, min_colour):
   # this is the basic function to get a colour given the range and domain
   # you can add a wrapper to shorten the parameter's length
   ratio = (this_IC - min_IC) / (max_IC - min_IC)
   result = [0,0,0]
   for i in range(3):
        result[i] = (max_colour[i] - min_colour[i]) * ratio + min_colour[i]
   return result

def translate_colour(colour):
    # from [255,255,255] to #FFFFFF
    result = '#'
    for c in colour:
        temp = hex(int(c))[2:]
        if len(temp) == 1:
            result += '0' + temp
        else:
            result += hex(int(c))[2:]
    return result

def IC(freq):
    # calculate information content
    return -math.log(freq)

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db=conn['hpo']
db = conn['uclex']
patient_db=conn['patients']
parser.add_option("--chrom",
                  dest="chrom",
                  help="which chrom to process?")
(options, args) = parser.parse_args()
version = 1
release = '2016_Aug'
CADD_cutoff = 15
description = 'merge all data and stats into one file. use Cadd_phred 15 as cutoff for p value calculation.'
genes = get_chrom_genes([options.chrom], db)
impute_file = '/cluster/project8/vyp/doug/uclex/uclex_phased.bed'
# get vcf file
vcf_location = '/SAN/vyplab/UCLex/current'
vcf_release = 'July2016'
vcf_file = os.path.join(vcf_location,'mainset_'+vcf_release+'_chr'+options.chrom+'.vcf.gz')
tbx = pysam.TabixFile(vcf_file) 
# get tbx header
tbx_header = []
for h in tbx.header:
    if h[:2] == '##': continue
    tbx_header = h.split('\t')
    break
patients_hpo_file = 'patients_hpo_snapshot_2016_Aug.tsv'
patients_hpo = parse_patients_hpo(patients_hpo_file)
hpo_freq_file = 'hpo_freq_'+release+'.json'
hpo_freq = json.load(open(hpo_freq_file,'r'))
related_pat = hpo_freq['related']['HP:0000001']
related_pat_a = len(related_pat)
unrelated_pat = hpo_freq['unrelated']['HP:0000001']
unrelated_pat_a = len(unrelated_pat)
# set up IC and colour
min_IC = 0
min_colour = [255, 0, 0]
max_colour = [198, 222, 255]
related_max_IC = IC(1./related_pat_a)
unrelated_max_IC = IC(1./unrelated_pat_a)
print 'related pat_a: %s, unrelated pat_a: %s' % (related_pat_a, unrelated_pat_a)
i = 0
#v = '7-138601865-G-A'
#this = vcf_query(variant_str=v)
#print this['het_samples']
#sys.exit()
for gene_id in genes:
    i += 1
    
    if not gene_id.startswith('ENSG'): gene_id = get_gene_by_name(db, gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    # already done?
    filename = os.path.join('gene_hpo',release,str(version),gene_id + '.json')
    inf = open(filename,'r')
    data = json.load(inf)
    inf.close()

    print '====='
    print gene_id
    print '%s / %s' % (i, len(genes))
    #print json.dumps(result, indent=4)
    #sys.exit()
    for mode in ['hom_comp','het']:
        # populate fisher p
        populate_fisher_p(data[mode],mode)
    outf = open(filename,'w')
    json.dump(data,outf,indent=4)
