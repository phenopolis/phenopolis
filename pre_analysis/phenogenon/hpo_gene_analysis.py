#!/bin/env python
'''
this script uses the results calculated by gene_hpo_analysis.py to rank the genes according to the p / phi values
'''
import json
import os
import sys
import time

release = '2016_Aug'
version = 4
hpo_freq_file = 'hpo_freq_'+release+'.json'
hpo_freq = json.load(open(hpo_freq_file,'r'))

path_to_files = os.path.join('gene_hpo',release,str(version))
result = {}#{hpo:{'release':release,'version':version,'cadd_cutoff':CADD_cutoff,'exac_recessive_cutoff':0.01,'exac_dominant_cutoff':0.001,'hpo_id':HP:0001,'data':{'related':{'dominant':[],'recessive':[]},'unrelated':{'dominant':[],'recessive':[]}}}}
for file in os.listdir(path_to_files):
    if file.endswith('.json'):
        print file
        file = os.path.join(path_to_files,file)
        data = json.load(open(file,'r'))
        for mode in ['hom_comp','het']:
            for hpo,value in data[mode].iteritems():
                result[hpo] = result.get(hpo,{
                    'hpo_id':hpo,
                    'cadd_cutoff':CADD_cutoff,
                    'exac_recessive_cutoff':0.01,
                    'exac_dominant_cutoff':0.001,
                    'release':release,
                    'version':version,
                    'data':{
                        'related':{
                            'dominant':[],
                            'recessive':[]
                            },
                        'unrelated':{
                            'dominant':[],
                            'recessive':[]
                        }}})
                for r in ['related','unrelated']:
                    if mode == 'hom_comp' and value[r+'_recessive_p_val'] <= 0.05:
                        result[hpo]['data'][r]['recessive'].append({
                            'gene_id':data['gene_id'],
                            'p_val':value[r+'_recessive_p_val'],
                        })
                    elif mode == 'het' and value[r+'_dominant_single_p_val'] <= 0.05:
                        result[hpo]['data'][r]['dominant'].append({
                            'gene_id':data['gene_id'],
                            'p_val':value[r+'_dominant_single_p_val'],
                        })

# sort based on p value
for hpo, value in result.iteritems():
    for r in ['related','unrelated']:
        for mode in ['dominant','recessive']:
            value['data'][r][mode] = sorted(value['data'][r][mode], key=lambda k: k['p_val'])

# write to file, each hpo stands for a file
for hpo, value in result.iteritems():
    outfile = os.path.join('hpo_gene',release,str(version),hpo+'.json')
    json.dump(value,open(outfile,'w'),indent=4)
