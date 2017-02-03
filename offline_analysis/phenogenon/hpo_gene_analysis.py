#!/bin/env python
'''
this script uses the results calculated by gene_hpo_analysis.py to rank the genes according to the p / phi values
'''
from __future__ import print_function, division
import json
import os
import sys
sys.path.append('../commons')
from phenopolis_utils import *
import time
from gene_hpo_analysis import *

# the gene-hpo stats are wrong. Need to re-calculate
release = '2016_Aug'
version = 7
hpo_freq_file = 'hpo_freq_'+release+'.json'
hpo_freq = json.load(open(hpo_freq_file,'r'))
result_file = ''
vp_ratio_cutoffs = {
        'r': float(OFFLINE_CONFIG['filters']['vp_ratio_rec']),
        'd': float(OFFLINE_CONFIG['filters']['vp_ratio_dom']),
        }
calc_cutoffs = {
    'exac_rec':0.01,
    'exac_hom':2,
    'exac_dom':0.0001,
    'cadd':15
}
path_to_files = os.path.join('gene_hpo',release,str(version))
# CADD cutoff 15 is default, and precalculated. Will rank for both related and unrelated
result = {}#{hpo:{'release':release,'version':version,'cadd_cutoff':CADD_cutoff,'exac_recessive_cutoff':0.01,'exac_dominant_cutoff':0.001,'hpo_id':HP:0001,'data':{'related':{'dominant':[],'recessive':[]},'unrelated':{'dominant':[],'recessive':[]}}}}

for file in os.listdir(path_to_files):
    if file.endswith('.json'):
        print(file)
        file = os.path.join(path_to_files,file)
        data = json.load(open(file,'r'))

        for mode in ['r','d']:
            v_p = populate_fisher_p(data, mode, calc_cutoffs)
            if not v_p['patients']: continue
            vp_ratio = len(v_p['variants']) / len(v_p['patients'])
            if vp_ratio < vp_ratio_cutoffs[mode]:continue
            for hpo,value in data['data'].items():
                if hpo == 'HP:0000001': continue
                result[hpo] = result.get(hpo,{
                    'hpo_id':hpo,
                    'cadd_cutoff':calc_cutoffs['cadd'],
                    'exac_recessive_cutoff':calc_cutoffs['exac_hom'],
                    'exac_dominant_cutoff':calc_cutoffs['exac_dom'],
                    'release':release,
                    'version':version,
                    'data':{
                        'unrelated':{
                            'dominant':[],
                            'recessive':[]
                        }}})
                if value[mode+'_unrelated_pass_p_val'] <= 0.05:
                    m = 'recessive' if mode == 'r' else 'dominant'
                    result[hpo]['data']['unrelated'][m].append({
                        'gene_id':data['gene_id'],
                        'p_val':value[mode+'_unrelated_pass_p_val'],
                    })

# sort based on p value
for hpo, value in result.items():
    for mode in ['dominant','recessive']:
        value['data']['unrelated'][mode] = sorted(value['data']['unrelated'][mode], key=lambda k: k['p_val'])

# write to file, each hpo stands for a file
for hpo, value in result.items():
    outfile_path = os.path.join('hpo_gene',release,str(version))
    mkdir_p(outfile_path)
    outfile = os.path.join(outfile_path,hpo+'.json')
    json.dump(value,open(outfile,'w'),indent=4)
