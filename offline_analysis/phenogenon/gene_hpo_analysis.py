#!/bin/env python
'''
Use KING to get independent individuals
/slms/gee/research/vyplab/UCLex/KING/KING/UCL-exome_unrelated.txt

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
import errno
#from plinkio import plinkfile

conn = pymongo.MongoClient(host='phenotips', port=27017)
db = conn['uclex']
hpo_db=conn['hpo']
patient_db=conn['patients']

debug = None

'''
defs
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def populate_fisher_p(data,mode,calc_cutoffs):
    # add fisher p value into the data.
    # also returns qualified patients / variants as in 
    # { patients: [], variants: []}
    # only counts unrelated patients and passed variants
    # if mode == het, has to calculate d_p_val for dominant p
    if not data: return None
    dubious_p = {}
    dom_single_p = {}
    dubious_single_p = {}
    qualified_p = {}
    qualified_v = []
    # get qualified patients
    for k1,v1 in data['patients'].items():
        count = 0 # controls qualified var count
        subcount = 0 # controls dubious var count
        for var in v1['variants']:
            if ((mode == 'r' and data['variants'][var]['exac_af'] <= calc_cutoffs['exac_rec'] and data['variants'][var]['exac_hom'] <= calc_cutoffs['exac_hom']) or (mode=='d' and data['variants'][var]['exac_af'] <= calc_cutoffs['exac_dom'])) and (not data['variants'][var]['cadd'] or data['variants'][var]['cadd'] >= calc_cutoffs['cadd']): 
                count += 1
                if 'filter' in data['variants'][var] and data['variants'][var]['filter'] != 'PASS':
                    subcount += 1
                else:
                    qualified_v.append(var)
        if (mode == 'r' and count > 1) or (mode == 'd' and count):
            qualified_p[k1] = {'unrelated':v1['unrelated']}
            if mode == 'd' and count == 1:
                dom_single_p[k1] = {'unrelated':v1['unrelated']}
            if subcount and ((mode == 'd' and count == subcount) or (mode == 'r' and count - subcount < 2)):
                dubious_p[k1] = {'unrelated':v1['unrelated']}
                dubious_single_p[k1] = {'unrelated':v1['unrelated']}
            if subcount and (mode == 'd' and count - subcount <= 1):
                dubious_single_p[k1] = {'unrelated':v1['unrelated']}
    # calculate pat_a
    mt = {'hom_comp':'recessive','het':'dominant'}
    pat_a = {
            'unrelated':data['unrelated_pat_a'],
            'related':data['related_pat_a'],
            'unrelated_single':data['unrelated_pat_a'],
            'related_single':data['related_pat_a'],
            'unrelated_pass':data['unrelated_pat_a'] - len([i for i,j in dubious_p.iteritems() if j['unrelated']]),
            'related_pass':data['related_pat_a'] - len(dubious_p),
            'unrelated_pass_single':data['unrelated_pat_a'] - len([i for i,j in dubious_single_p.iteritems() if j['unrelated']]),
            'related_pass_single':data['related_pat_a'] - len(dubious_single_p)
            }
    # pat_g
    related_pat_g = len(qualified_p)
    unrelated_pat_g = len([i for i,j in qualified_p.items() if j['unrelated']])
    unrelated_pass  = set([i for i,j in qualified_p.items() if j['unrelated']]) - set(dubious_p.keys())
    pat_g = {
            'related':related_pat_g,
            'unrelated':unrelated_pat_g,
            'related_pass':related_pat_g - len(dubious_p),
            'unrelated_pass':unrelated_pat_g - (pat_a['unrelated'] - pat_a['unrelated_pass']),
            'unrelated_pass':len(unrelated_pass),
            'related_single':len(dom_single_p),
            'unrelated_single':len([i for i,j in dom_single_p.items() if j['unrelated']]),
            'related_pass_single':len(set(dom_single_p.keys()) - set(dubious_single_p.keys())),
            'unrelated_pass_single':len(set([i for i,j in dom_single_p.iteritems() if j['unrelated']]) - set([i for i,j in dubious_single_p.items() if j['unrelated']]))
            }
    # populate result
    result = {
            'patients':list(unrelated_pass),
            'variants':{},
            }
    for p in unrelated_pass:
        result['variants'].update({v:1 for v in data['patients'][p]['variants'] if v in qualified_v})
    # calculate
    for k1,v1 in data['data'].items():
        if k1 == 'HP:0000001':continue
        # get pat_h
        unrelated_pat_h = v1['unrelated_pat_h']
        related_pat_h = v1['related_pat_h']
        pat_h = {
                'unrelated':unrelated_pat_h,
                'related':related_pat_h,
                'unrelated_single':unrelated_pat_h,
                'related_single':related_pat_h,
                'unrelated_pass':unrelated_pat_h - len(set([i for i in v1['p'][mode] if data['patients'][i]['unrelated']]) & set(dubious_p.keys())),
                'related_pass':related_pat_h - len(set(v1['p'][mode]) & set(dubious_p.keys())),
                'unrelated_pass_single':unrelated_pat_h - len(set([i for i in v1['p']['d'] if data['patients'][i]['unrelated']]) & set(dubious_single_p.keys())),
                'related_pass_single':related_pat_h - len(set(v1['p']['d']) & set(dubious_single_p.keys()))
                }
        # get pat_gh
        pat_gh = {
                'unrelated':len(set([i for i in v1['p'][mode] if data['patients'][i]['unrelated']]) & set(qualified_p.keys())),
                'related':len(set(v1['p'][mode]) & set(qualified_p.keys())),
                'unrelated_pass':len((set([i for i in v1['p'][mode] if data['patients'][i]['unrelated']]) & set(qualified_p.keys())) - set(dubious_p.keys())),
                'related_pass':len((set(v1['p'][mode]) & set(qualified_p.keys())) - set(dubious_p.keys())),
                'unrelated_single':len(set([i for i in v1['p']['d'] if data['patients'][i]['unrelated']]) & set(dom_single_p.keys())),
                'related_single':len(set(v1['p']['d']) & set(dom_single_p.keys())),
                'unrelated_pass_single':len((set([i for i in v1['p']['d'] if data['patients'][i]['unrelated']]) & set(dom_single_p.keys())) - set(dubious_single_p.keys())),
                'related_pass_single':len((set(v1['p']['d']) & set(dom_single_p.keys())) - set(dubious_single_p.keys()))
                }
        #calculate
        if mode == 'r':
            # for key in ['related','unrelated','related_pass','unrelated_pass']:
            # only use unrelated_pass, rest can be calculated on live using js
            for key in ['unrelated_pass']:
                v1['r_'+key+'_p_val'] = fisher.pvalue(pat_a[key]-pat_h[key]-pat_g[key]+pat_gh[key],pat_h[key]-pat_gh[key],pat_g[key]-pat_gh[key],pat_gh[key]).right_tail
        else:
            if debug and k1 == debug:
                tk ='unrelated'
                print pat_a[tk]
                print pat_g[tk]
                print pat_h[tk]
                print pat_gh[tk]
            #for key in pat_gh:
            for key in ['unrelated_pass','unrelated_pass_single']:
                v1['d_'+key+'_p_val'] = fisher.pvalue(pat_a[key]-pat_h[key]-pat_g[key]+pat_gh[key],pat_h[key]-pat_gh[key],pat_g[key]-pat_gh[key],pat_gh[key]).right_tail

    return result

def get_gene_hpo(gene_id):
    # similar to get_rare_var_p_hpo, but this one has hpo as key 
    #return {'hom_comp':{'HP:0001234':{id:'HP:0001234', name:'hell', data:{p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}}},
    #        'het':{'HP:0001234':{id:'HP:0001234', name:'hell', data:{p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}}}}

    # sometimes variant is not in vcf. move it to debug/bad_variants for inspection and later clean
    # for het, cut at exac_af = 0.001
    # get all variants on this gene
    #all_vars = db.genes.find_one({'gene_id':gene_id})['variant_ids']
    # when all_vars too long, cursor will die. convert to list
    all_vars = list(db.variants.find({'genes':gene_id}))
    results = {'data':{}, 'patients':{}, 'variants':{}} 
    #for v in all_vars:
    for var in all_vars:
        v = var['variant_id']
        '''
        v=None
        if 'variant_id' not in var:
            # some var has 'VARIANT_ID' instead of 'variant_id', correct it
            v=var['VARIANT_ID']
            db.variants.update({'VARIANT_ID':v},{'$rename':{'VARIANT_ID':'variant_id'}})
        else:
            v = var['variant_id']
        '''
        # get most_severe_consequence, remove synonymous/intronic ones
        # refer to http://www.ensembl.org/info/genome/variation/predicted_data.html
        # if most severe impact is 'MODIFIER' or ('LOW' and not 'splice'), continue
        low = 1
        for con in [cons for cons in var['transcript_consequences'] if cons['gene_id'] == gene_id]:
            modifier = con.get('impact','MODIFIER')
            if modifier not in ['MODIFIER','LOW']:
                low = 0
                break
            elif modifier == 'LOW':
                if 'splice_region_variant' in con.get('consequence_terms',[]):
                    low = 0
                    break
        if low: continue
        try:
            var_obj = orm.Variant(variant_id=v,db=db)
        except:
            print v+'not in vcf'
            continue
        # failed the filter?
        if var_obj.filter=='FAIL': 
            continue
        var['cadd'] = max(var_obj.cadd) if type(var_obj.cadd) is list else var_obj.cadd
        cadd_phred = var['cadd']

        # get exac_af/hom
        exac_af = exac_hom = 0
        if var_obj.in_exac:
            # exac_af = var_obj.EXAC['AF'], gives error on multiallelic site, such as 8-43002117-C-A
            if var_obj.ExAC_freq:
                exac_af = 0. if not var_obj.ExAC_freq['total_ans'] else float(var_obj.ExAC_freq['total_acs'])/var_obj.ExAC_freq['total_ans']
                exac_hom = var_obj.ExAC_freq.get('total_homs',0)
        # not interested if af is higher than exac_recessive cutoff
        if exac_af > filter_cutoffs['exac_rec'] or exac_hom > filter_cutoffs['exac_hom']:
            continue
        '''
        var_obj.ExAC_freq['total_homs'] encode the number of homs. good for recessive cases. will incorporate
        '''
        # get relevant info from vcf
        #uclex_af = var_obj.allele_freq # not using it
        
        # dealing with hom patients. also add it to het.
        # will need to deal with both_het !!! their af are different!!!
        hom_p = var_obj.hom_samples
        for p in hom_p:
            if p not in patients_hpo['related']: continue
            populate_mode_p(results, 'r', p, exac_af, exac_hom, v, cadd_phred, var_obj.filter)

        # dealing with het patients. note to check length of exac_af. longer than one?
        # also added it to 'hom_comp'
        het = var_obj.het_samples
        for p in het:
            if p not in patients_hpo['related']: continue
            populate_mode_p(results, 'd', p, exac_af, exac_hom, v, cadd_phred, var_obj.filter)
    # prepare to return the result
    results['version'] = version
    results['description'] = ''
    results['gene_id'] = gene_id
    results['release'] = release
    results['related_pat_a'] = related_pat_a
    results['unrelated_pat_a'] = unrelated_pat_a
    return results

def populate_mode_p(results, mode, p, exac_af, exac_hom, v, cadd_phred, filter):
    # get all hpos of the patient. Note that related has all the patients in there
    hpos = patients_hpo['related'][p]['hpo']
    minimum_set = hpo_minimum_set(hpo_db, hpo_ids=hpos)
    minimum_set_array = list(hpo_db.hpo.find({'id':{'$in':minimum_set}}))
    hpos_names = [i['name'][0] for i in minimum_set_array]
    unrelated = 0
    if p in patients_hpo['unrelated']:
        unrelated = 1
    # add patients
    if p not in results['patients']:
        results['patients'][p] = {'hpo':hpos_names, 'unrelated':unrelated, 'variants':[]}
    # add variants according to homozygosity
    # hom = 2, het = 1
    copy = 2 if mode == 'r' else 1
    results['patients'][p]['variants'].extend([v]*copy)
    if v not in results['variants']:
        results['variants'][v] = {
                'exac_af':exac_af,
                'filter':filter,
                'exac_hom':exac_hom,
                'cadd':cadd_phred
        }
    
    for h in hpos:
        hpo_dict = hpo_db.hpo.find_one({'id':h})
        if h not in results['data']:
            # initialise
            results['data'][h] = {
                    'id':h,
                    'name':hpo_dict['name'][0],
                    'is_a':hpo_dict.get('is_a',[]),
                    'unrelated_pat_h':len(hpo_freq['unrelated'].get(h,[])),
                    'related_pat_h':len(hpo_freq['related'][h]),
                    'p':{'d':[],'r':[]}
                    }
        if p in results['data'][h]['p'][mode]:
            if mode == 'd' and p not in results['data'][h]['p']['r']:
                # dominant mode and more than 1 variants, copy it to recessive
                results['data'][h]['p']['r'].append(p)
        else:
            results['data'][h]['p'][mode].append(p)
        if mode == 'r' and exac_af <= filter_cutoffs['exac_dom'] and p not in results['data'][h]['p']['d']:
            results['data'][h]['p']['d'].append(p)

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
'''
main
'''
if __name__ == "__main__":
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--chrom",
                      dest="chrom",
                      help="which chrom to process?")
    (options, args) = parser.parse_args()
    version = 6
    release = '2016_Aug'
    # filter_cutoffs is used to include patients and variants in the result file
    # calc_cutoffs is used to get patients for phenogenon test
    filter_cutoffs = {
        'exac_rec':0.01,
        'exac_hom':10,
        'exac_dom':0.001,
        'cadd':0
    }
    calc_cutoffs = {
        'exac_rec':0.01,
        'exac_hom':2,
        'exac_dom':0.0001,
        'cadd':15
    }

    description = 'add missing filters\ndo not keep cleaned variant_id\nremove redundant patient/var info'
    genes = get_chrom_genes([options.chrom], db)
    #genes=['TRIM50'];
    #impute_file = '/cluster/project8/vyp/doug/uclex/uclex_phased.bed'
    patients_hpo_file = 'patients_hpo_snapshot_'+release+'.tsv'
    patients_hpo = parse_patients_hpo(patients_hpo_file)
    hpo_freq_file = 'hpo_freq_'+release+'.json'
    hpo_freq = json.load(open(hpo_freq_file,'r'))
    related_pat = hpo_freq['related']['HP:0000001']
    related_pat_a = len(related_pat)
    unrelated_pat = hpo_freq['unrelated']['HP:0000001']
    unrelated_pat_a = len(unrelated_pat)
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
        filepath = os.path.join('gene_hpo', release, str(version))
        mkdir_p(filepath)
        filename = os.path.join(filepath, gene_id + '.json')
        if os.path.isfile(filename) and os.stat(filename).st_size:
            bad = 0 # change it to 1 if want to overwrite anyway.
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
        result = get_gene_hpo(gene_id)
        #print json.dumps(result, indent=4)
        #sys.exit()
        #for mode in ['r','d']:
            # populate fisher p
            # populate_fisher_p(result,mode,calc_cutoffs)
        outf = open(filename,'w')
        json.dump(result,outf,indent=4)


            # if you want a static dot file, run the following
            #dot_data = transform(rare_p_hpo[mode], hpo_freq)
            #test = draw_dot_graph(dot_data)
            #outf.write(test)
