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
from scipy.stats import binom
import os
import sys
import time
import rest

conn = pymongo.MongoClient(host='localhost', port=27017)
hpo_db=conn['hpo']
db = conn['uclex-old']
patient_db=conn['patients']
###############################
# parse options
###############################
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-m", "--mode", dest="mode", default='nearest',
                  help='either nearest or average [default: %default]')
                  
(options, args) = parser.parse_args()

# vars that don't have proper consequence terms will be thrown here. most of them come from X chrom
badvars = open('bad_vars.txt','w')
'''
defs
'''

def get_rare_var_hpo_p(gene_id, db, patient_db, hpo_freq, unrelated, version):
    # similar to get_rare_var_p_hpo, but this one has hpo as key 
    #return {'hom_comp':{'HP:0001234':{id:'HP:0001234', name:'hell', data:{p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}}},
    #        'het':{'HP:0001234':{id:'HP:0001234', name:'hell', data:{p_id1:[{v_id,exac_af,cadd_phred}],p_id2:[{v_id,exac_af,cadd_phred}]}}}}

    # sometimes variant is not in vcf. move it to debug/bad_variants for inspection and later clean
    # for het, cut at exac_af = 0.001
    bad_var_file = open('views/debug/bad_variants', 'w')
    # get all variants on this gene
    #all_vars = db.genes.find_one({'gene_id':gene_id})['variant_ids']
    # when all_vars too long, cursor will die. convert to list
    all_vars = list(db.variants.find({'genes':gene_id}))
    results = {'hom_comp':{}, 'het':{}} 
    #for v in all_vars:
    for var in all_vars:
        v=None
        if 'variant_id' not in var:
            # some var has 'VARIANT_ID' instead of 'variant_id', correct it
            v=var['VARIANT_ID']
            db.variants.update({'VARIANT_ID':v},{'$rename':{'VARIANT_ID':'variant_id'}})
        else:
            v = var['variant_id']

        # get consequence, remove synonymous ones
        if 'Consequence' not in var:
            #borrow vep.Consequence
            if 'vep_annotations' not in var:
                print v
                badvars.write(v+'\n')
                continue
                chrom,pos,ref,alt = v.split('-')
                returned = rest.vep_anno(chrom,pos,ref,alt)
                #print returned[0]
            var['Consequence'] = [t['Consequence'] for t in var['vep_annotations']]
            db.variants.update({'variant_id':v},{'$set':{'Consequence':var['Consequence']}})
        if not set(var['Consequence']) - set(['synonymous_variant','upstream_gene_variant','downstream_gene_variant','regulatory_region_variant','MODIFIER','intron_variant','non_coding_transcript_variant','intron_variant&non_coding_transcript_variant','non_coding_transcript_exon_variant&non_coding_transcript_variant']):
            continue
        # get cadd phred scores. if not exist, assign 50
        try:
            var_obj = orm.Variant(variant_id=v,db=db)
        except:
            continue
        #print var_obj.consequence
        if 'cadd' not in var:
            for attempt in range(10):
                try:
                    # sometimes the url connection is bad
                    var['cadd'] = var_obj.cadd
                    break
                except:
                    time.sleep(5)
        if 'cadd' not in var:
            print 'tried too many times, abort connecting to url'
            raise
        if isinstance(var['cadd'], float) or isinstance(var['cadd'], int):
            cadd_phred = var['cadd']
        else:
            cadd_phred = None if not var['cadd'] else var['cadd'].get('phred',None)

        # get exac_af
        exac_af = 0
        if 'in_exac' not in var:
            rest.exac_anno(v)
            var = db.variants.find_one({'variant_id':v})
            #if 'allele_freq' not in VAR:
            #    print '====allele_freq not in external anno==='
            #    print VAR
            #elif VAR['allele_freq']:
            #    exac_af = float(VAR['allele_freq'])
        #elif var['in_exac']:
        if var['in_exac']:
            #if 'allele_freq' not in var['EXAC']:
            #    VAR = rest.exac_anno(v)
            #    if 'allele_freq' not in VAR:
            #        print '====allele_freq not in external anno==='
            #        print VAR
            #    elif VAR['allele_freq']:
            #        exac_af = float(VAR['allele_freq'])
            #else:
            if not var['EXAC']['allele_freq']:
                exac_af = 0.
            else:
                exac_af = float(var['EXAC']['allele_freq'])
        # not interested if af is > 0.01
        if exac_af > 0.01:
            continue
        # get relevant info from vcf
        uclex_af = var_obj.allele_freq # not using it

        # dealing with hom patients. also add it to het.
        # will need to deal with both_het !!! their af are different!!!
        hom_p = var_obj.hom_samples
        for p in hom_p:
            #if p not in unrelated:
            #    continue
            populate_mode_p(results, ['hom_comp', 'het'], p, exac_af, patient_db, hpo_freq, v, cadd_phred)

        # dealing with het patients. note to check length of exac_af. longer than one?
        # also added it to 'hom_comp'
        het = var_obj.het_samples
        for p in het:
            #if p not in unrelated:
            #    continue
            populate_mode_p(results, ['het'], p, exac_af, patient_db, hpo_freq, v, cadd_phred)
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
    return {'hom_comp':{'data':results['hom_comp'],'version':version},'het':{'data':results['het'],'version':version}}


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

def populate_mode_p(results, modes, p, exac_af, patient_db, hpo_freq, v, cadd_phred):
    hpos = get_patient_observed_hpo(p, patient_db)
    total = hpo_freq['HP:0000001']['size']
    max_IC = IC(1./total)
    min_IC = 0
    min_colour = [255, 0, 0]
    max_colour = [198, 222, 255]
    if hpos:
        all_hpos = []
        # get all its ancestors
        for h in hpos:
            if not h[0]:
                continue
            all_hpos += get_hpo_ancestors_array(hpo_db, h[0])
        all_hpos = set(all_hpos)
        for mode in modes:
            # copy if mode = hom_comp
            copy = 2 if mode == 'hom_comp' else 1
            for h in all_hpos:
                hpo_dict = hpo_db.hpo.find_one({'id':h})
                if h not in results[mode]:
                    # initialise
                    this_IC = IC(hpo_freq[h]['freq'])
                    colour = get_gradient_colour_raw(this_IC, max_IC, min_IC, max_colour, min_colour)
                    results[mode][h] = {'id':h, 'name':hpo_dict['name'][0], 'is_a':hpo_dict.get('is_a',[]), 'fillcolor':translate_colour(colour), 'observed_freq':None, 'expect_freq':hpo_freq[h]['raw'],  'data':{}}
                if p in results[mode][h]['data']:
                    results[mode][h]['data'][p]['var'].extend([{'variant_id':v,'exac':exac_af, 'cadd':cadd_phred}]*copy)
                    if len(modes) == 1:
                        # het mode, and this patient has more than one rare variants
                        # copy it to hom_comp
                        this_IC = IC(hpo_freq[h]['freq'])
                        colour = get_gradient_colour_raw(this_IC, max_IC, min_IC, max_colour, min_colour)
                        results['hom_comp'][h] = results['hom_comp'].get(h,{'id':h, 'name':hpo_dict['name'][0], 'is_a':hpo_dict.get('is_a',[]), 'fillcolor':translate_colour(colour),'observed_freq':None, 'expect_freq':hpo_freq[h]['raw'],  'data':{}})

                        if p in results['hom_comp'][h]['data']:
                            results['hom_comp'][h]['data'][p]['var'].append({'variant_id':v,'exac':exac_af, 'cadd':cadd_phred})
                        else:
                            results['hom_comp'][h]['data'][p] = results['het'][h]['data'][p]
                else:
                    results[mode][h]['data'][p] = {'hpo':[t[1] for t in hpos], 'var':[{'variant_id':v, 'exac':exac_af, 'cadd':cadd_phred}]*copy}

def get_chrom_genes(chroms, db):
    # give chrom numbers, get all genes on them
    result = []
    for chrom in chroms:
        genes = [g['gene_id'] for g in db.genes.find({'chrom':str(chrom)})]
        result.extend(genes)
    return result

'''
main
'''
version = 2
genes = get_chrom_genes([2], db)
#genes=[g['gene_id'] for g in db.genes.find(fields={'_id':False})]
hpo_freq = get_hpo_size_freq('hpo_freqi.tsv')
unrelated = open('/slms/gee/research/vyplab/UCLex/mainset_July2016/kinship/UCL-exome_unrelated.txt','r').readlines()
unrelated = [i.rstrip() for i in unrelated]

#genes=['GUCY2D'];
i = 0
#v = '7-138601865-G-A'
#this = vcf_query(variant_str=v)
#print this['het_samples']
#sys.exit()
for gene_id in genes:
    i += 1
    #if gene_id in ['ENSG00000213930','ENSG00000258728','ENSG00000137070','ENSG00000188536','ENSG00000225323','ENSG00000206172']:
        # having trouble with this gene. it is on chr9, but somehow db.variants chr22 vars are linked to this gene, such as 22-17060797-T-C
    #    continue
    if not gene_id.startswith('ENSG'): gene_id = get_gene_by_name(db, gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    # already done?
    hom_file = os.path.join('static','dot',gene_id + '_hom_comp.json')
    het_file = os.path.join('static','dot',gene_id + '_het.json')
    if os.path.isfile(hom_file) and os.stat(hom_file).st_size and os.path.isfile(het_file) and os.stat(het_file).st_size:
        hom_json = json.load(open(hom_file,'r'))
        het_json = json.load(open(het_file,'r'))
        if hom_json.get('version', None) == version and het_json.get('version', None) == version:
            continue

    print '====='
    print gene_id
    print '%s / %s' % (i, len(genes))
    rare_p_hpo = get_rare_var_hpo_p(gene_id, db, patient_db, hpo_freq, unrelated, version)
    for mode in ['hom_comp','het']:
        filename = os.path.join('static','dot',gene_id + '_' + mode + '.json')
        #filename = os.path.join(gene_id + '_' + mode + '.json')
        outf = open(filename,'w')
        json.dump(rare_p_hpo[mode],outf)

        # if you want a static dot file, run the following
        #dot_data = transform(rare_p_hpo[mode], hpo_freq)
        #test = draw_dot_graph(dot_data)
        #outf.write(test)
