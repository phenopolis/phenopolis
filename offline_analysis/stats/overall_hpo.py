#!/bin/env python
import sys
import pymongo
import math
import os
import json
import re
sys.path.append('..')
#from lookups import *

'''
generagte hpo graph for all valid patients in uclex
{related:{
    HP:00001:{
        count:1000,
        freq:0.4, 
        name:foo,
        is_a:HP:0001,
        cohort:{
            UKIRDC:{
                M:{
                    count:100,
                    solved:10,
                    candidate:20
                },
                F:{},
                U:{},
            },
            WEBSTER:{}...
        }
    },
unrelated:{}
}
one can then extract the overall stats from HP:0000001
'''

conn = pymongo.MongoClient(host='phenotips', port=27017)
hpo_db = conn['hpo']
db = conn['uclex']
patient_db = conn['patients']


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
version = 2
release = '2016_Aug'
patients_hpo_file = '../patients_hpo_snapshot_'+release+'.tsv'
hpo_freq_file = '../hpo_freq_'+release+'.json'
outf = open('overall_hpo_'+release+'_'+str(version)+'.json','w')
hpo_freq = json.load(open(hpo_freq_file,'r'))

# get the hpo snapshot
patients_hpo = parse_patients_hpo(patients_hpo_file)

# only getting related patients will be enough. since they include unrelated patients
'''
{p_id:{cohort:UKIRDC,sex:[M,F,U],solved=[1=candidate,2=solved,0=unsolved]}...}
'''
    
# get patients details
patients = {}
for p in patients_hpo['related']:
    temp = p
    # correcting some patient id
    if p == 'Vulliamy_UCLG.177': temp = 'Vulliamy_UCLG-177'
    if p[:14] == 'Vulliamy_UCLG.':
        temp = '-'.join(p.split('.'))
    elif p[:16] == 'Vulliamy_May2015' and len(p.split('_')) == 3:
        fields = p.split('_')
        temp = '_'.join(fields[:2]+fields[:2]+[fields[2]])
    elif p[:18] == 'Vulliamy_April2016' and re.match(r'\d',p.split('_')[2]):
        fields = p.split('_')
        temp = '_'.join(fields[:2]+['AW'+fields[2]])
    elif p[:22] == 'Vulliamy_April2016_AW_':
        fields = p.split('_')
        temp = '_'.join(fields[:3]) + fields[3]
    pat = patient_db.patients.find_one({'external_id':temp})
    if not pat:
        pat = patient_db.patients.find_one({'external_id':re.compile(temp)})
    solved = 2
    print p
    if pat['solved']['status']=='unsolved':
        solved = 0
        if pat.get('genes',[]):
            solved = 1
    # Jing Yu and Nikolas Pontikos are actually UKIRDC
    cohort = pat['contact']['name']
    if cohort in ['Jing Yu']:
        cohort = 'UKIRDC'
    
    # get sex. there are some patients having there sex coded as 'O'...female??
    sex = pat['sex']
    if sex == 'O': sex = 'U'
    patients[p] = {
            'cohort':cohort,
            'sex':sex,
            'solved':solved
            }
# populate result
result = {'related':{},'unrelated':{}}

for mode in ['related','unrelated']:
    total = len(hpo_freq[mode]['HP:0000001'])
    for hpo,pats in hpo_freq[mode].iteritems():
        hpo_dict = hpo_db.hpo.find_one({'id':hpo})
        count = len(pats)
        result[mode][hpo] = {
            'count':count,
            'freq':float(count)/total,
            'name':hpo_dict['name'][0],
            'is_a':hpo_dict.get('is_a',[]),
            'cohort':{}
            }
        for p in pats:
            cohort = patients[p]['cohort']
            result[mode][hpo]['cohort'][cohort] = result[mode][hpo]['cohort'].get(cohort,{
                'M':{'unsolved':0,'solved':0,'candidate':0},
                'F':{'unsolved':0,'solved':0,'candidate':0},
                'U':{'unsolved':0,'solved':0,'candidate':0}
                })
            temp = result[mode][hpo]['cohort'][cohort][patients[p]['sex']]
            if patients[p]['solved'] == 1:
                temp['candidate'] += 1
            elif patients[p]['solved'] == 2:
                temp['solved'] += 1
            else:
                temp['unsolved'] += 1
    

# write to file
json.dump({'release':release,'version':version,'data':result},outf,indent=4)
