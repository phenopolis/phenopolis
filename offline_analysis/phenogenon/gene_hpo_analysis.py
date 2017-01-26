#!/bin/env python
'''
result = 'digraph {
6 [style="filled", fixedsize="true", fontsize="15", shape="circle", width="0.75", fillcolor="powderblue", label="Retina", data={p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}, expect_freq="23/2345", observed_freq="123/2345", id="HP:0000479", color="transparent"];
7 [style="filled", fixedsize="true", fontsize="15", shape="circle", width="0.75", fillcolor="powderblue", label="Retinal Dystrophy", data={p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}, expect_freq="23/2345", observed_freq="123/2345", id="HP:0000556", color="transparent"];
6 -> 7 [color="#000000", lty="solid"];
6 -> 39 [color="#000000", lty="solid"];
6 -> 105 [color="#000000", lty="solid"];
7 -> 112 [color="#000000", lty="solid"];
7 -> 114 [color="#000000", lty="solid"];
7 -> 172 [color="#000000", lty="solid"];
39 -> 172 [color="#000000", lty="solid"];
39 -> 183 [color="#000000", lty="solid"];
39 -> 40 [color="#000000", lty="solid"];
}

width's square propotional to observed/expected sample size. zero observed width = 0.01
filled color is propotional to IC of the HPO term, ranging from #FFFFFF to #C6DEFF
equivelent to from (255,255,255) to (198, 222, 255)

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
from plinkio import plinkfile

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db=conn['hpo']
db = conn['uclex']
patient_db=conn['patients']

# vars that don't have proper consequence terms will be thrown here. most of them come from X chrom
badvars = open('bad_vars.txt','w')
'''
defs
'''
def check_var(var, patient_array):
    # leave this for the time being!
    # first use Doug's result to impute the missed vars from available individuals, then use the existing ratios to impute the rest missing data. Note that some vars don't have imputation data
    # return hpos that are bad for this var
    missed_nohpo = called_nohpo = missed_hpo = called_hpo = 0
    chrom,pos,ref,alt = var.split('-')
    record = tbx.fetch(str(chrom),int(pos)-1,int(pos))
    missed_patients = []
    for row in record:
        row=row.split('\t')
        # genotype starts from row[9]
        for ind in range(9,len(tbx_header)):
            if tbx_header[ind] not in patient_array: continue
            called = 0 if row[ind][0]==row[2]=='.' else 1
            if not called:
                missed_patients.append(tbx_header[ind])
    print patient_array
    print var
    sys.exit()
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
    print obs
    if missed_nohpo == missed_hpo == 0: return 1
    # use chi square without Yates' correction to be conserved
    p_val = chi2_contingency(obs,correction=False)[1]
    return p_val

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
                if variant['exac'] <= exac_cutoff and ((not variant['cadd']) or variant['cadd'] >= CADD_cutoff): count += 1
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
                    if variant['exac'] <= exac_cutoff and ((not variant['cadd']) or variant['cadd'] >= CADD_cutoff): count += 1
                if count > 1:
                    failed.append(p)
            for p in failed:
                del v['data'][p]
        write_to_data(data,dp,'dominant_single')
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
    results = {'hom_comp':{}, 'het':{}} 
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

        # get exac_af
        exac_af = 0
        if var_obj.in_exac:
            # exac_af = var_obj.EXAC['AF'], gives error on multiallelic site, such as 8-43002117-C-A
            if var_obj.ExAC_freq:
                exac_af = 0. if not var_obj.ExAC_freq['total_ans'] else float(var_obj.ExAC_freq['total_acs'])/var_obj.ExAC_freq['total_ans']
        # not interested if af is > 0.01
        if exac_af > 0.01:
            continue
        # get relevant info from vcf
        #uclex_af = var_obj.allele_freq # not using it
        
        # dealing with hom patients. also add it to het.
        # will need to deal with both_het !!! their af are different!!!
        hom_p = var_obj.hom_samples
        for p in hom_p:
            if p not in patients_hpo['related']: continue
            populate_mode_p(results, ['het','hom_comp'], p, exac_af, v, cadd_phred, var_obj.filter)

        # dealing with het patients. note to check length of exac_af. longer than one?
        # also added it to 'hom_comp'
        het = var_obj.het_samples
        for p in het:
            if p not in patients_hpo['related']: continue
            populate_mode_p(results, ['het'], p, exac_af, v, cadd_phred, var_obj.filter)
    # clean het with cut at exac_af==0.001
    print 'before clean, rare_het length = %s' % len(results['het'])
    bad_h = []
    for key, h in results['het'].iteritems():
        bad = 1
        for p, exac in h['data'].iteritems():
            if min([e['exac'] for e in exac['var']]) <= 0.001:
                bad = 0
        if bad:
            bad_h.append(key)
    for h in bad_h:
        del results['het'][h]
    # add version
    #for mode in ['hom_comp','het']:
    #    results[mode]['version'] = version
    print 'after clean, rare_het length = %s' % len(results['het'])
    print 'rare_hom_comp length = %s' % len(results['hom_comp'])
    return {'version':version, 'description': '', 'gene_id':gene_id,'release':release,'hom_comp':results['hom_comp'],'het':results['het']}


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

def populate_mode_p(results, modes, p, exac_af, v, cadd_phred, filter):
    # get all hpos of the patient. Note that related has all the patients in there
    hpos = patients_hpo['related'][p]['hpo']
    minimum_set = hpo_minimum_set(hpo_db, hpo_ids=hpos)
    minimum_set_array = list(hpo_db.hpo.find({'id':{'$in':minimum_set}}))
    hpos_names = [i['name'][0] for i in minimum_set_array]
    unrelated = 0
    if p in patients_hpo['unrelated']:
        unrelated = 1

    for mode in modes:
        # copy if mode = homozygous
        copy = 2 if len(modes) == 2 else 1
        for h in hpos:
            hpo_dict = hpo_db.hpo.find_one({'id':h})
            this_related_IC = IC(len(hpo_freq['related'][h])/float(related_pat_a))
            related_colour = get_gradient_colour_raw(this_related_IC, related_max_IC, min_IC, max_colour, min_colour)
            this_unrelated_IC = IC(len(hpo_freq['unrelated'][h])/float(unrelated_pat_a)) if h in hpo_freq['unrelated'] else 0
            unrelated_colour = get_gradient_colour_raw(this_unrelated_IC, related_max_IC, min_IC, max_colour, min_colour)
            if h not in results[mode]:
                # initialise
                results[mode][h] = {
                        'id':h,
                        'name':hpo_dict['name'][0],
                        'is_a':hpo_dict.get('is_a',[]),
                        'related_colour':translate_colour(related_colour),
                        'unrelated_colour':translate_colour(unrelated_colour),
                        'unrelated_pat_a':unrelated_pat_a,
                        'related_pat_a':related_pat_a,
                        'unrelated_pat_h':len(hpo_freq['unrelated'].get(h,[])),
                        'related_pat_h':len(hpo_freq['related'][h]),
                        'data':{}
                        }
            
            # populate data
            if p in results[mode][h]['data']:
                results[mode][h]['data'][p]['var'].extend([{'variant_id':v,'exac':exac_af, 'cadd':cadd_phred}]*copy)
                if len(modes) == 1:
                    # het mode, and this patient has more than one rare variants
                    # copy it to hom_comp
                    results['hom_comp'][h] = results['hom_comp'].get(h,{
                        'id':h,
                        'name':hpo_dict['name'][0],
                        'is_a':hpo_dict.get('is_a',[]),
                        'related_colour':translate_colour(related_colour),
                        'unrelated_colour':translate_colour(unrelated_colour),
                        'unrelated_pat_a':unrelated_pat_a,
                        'related_pat_a':related_pat_a,
                        'unrelated_pat_h':len(hpo_freq['unrelated'].get(h,[])),
                        'related_pat_h':len(hpo_freq['related'][h]),
                        'data':{}}
                        )

                    if p in results['hom_comp'][h]['data']:
                        results['hom_comp'][h]['data'][p]['var'].append({'variant_id':v,'exac':exac_af, 'cadd':cadd_phred, 'filter':filter})
                    else:
                        results['hom_comp'][h]['data'][p] = results['het'][h]['data'][p]
            else:
                results[mode][h]['data'][p] = {'hpo':hpos_names, 'unrelated':unrelated, 'var':[{'variant_id':v, 'exac':exac_af, 'cadd':cadd_phred, 'filter':filter}]*copy}

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
    version = 2
    release = '2016_Aug'
    CADD_cutoff = 15
    description = 'version 1:merge all data and stats into one file;version 2:correct a bug when one variant is on multiple genes with messed up consequence annotation'
    genes = get_chrom_genes([options.chrom], db)
    #genes=['TRIM50'];
    #impute_file = '/cluster/project8/vyp/doug/uclex/uclex_phased.bed'
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
    patients_hpo_file = 'patients_hpo_snapshot_'+release+'.tsv'
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
        for mode in ['hom_comp','het']:
            # populate fisher p
            populate_fisher_p(result[mode],mode)
        outf = open(filename,'w')
        json.dump(result,outf,indent=4)

            # if you want a static dot file, run the following
            #dot_data = transform(rare_p_hpo[mode], hpo_freq)
            #test = draw_dot_graph(dot_data)
            #outf.write(test)
