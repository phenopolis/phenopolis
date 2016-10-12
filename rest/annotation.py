
"""
Variant and gene annotation.
Will use pymongo db as cache to avoid unecessary networks access.

Sources:
    VEP
    RVS
    ExAC
    myvariant.info
"""

import requests
import re
import json
import pymongo
# myvariant.info
import myvariant
# mygene.info
import mygene
# RVS
from varnorm.varcharkey import VarCharKey


mv = myvariant.MyVariantInfo()

mg = mygene.MyGeneInfo()

def ensembl_region(region,features='gene'):
    #'feature=gene;feature=transcript;feature=cds;feature=exon'
    r=requests.get('http://grch37.rest.ensembl.org/overlap/region/human/{}?feature={}'.format(region,features),headers={'Content-Type':'application/json'})
    return r.json()

def ensembl_xrefs(ID):
    r=requests.get('http://grch37.rest.ensembl.org/xrefs/id/{}'.format(ID),headers={'Content-Type':'application/json'})
    return r.json()

def mygene_anno(gene_symbol=None,gene_id=None, values='name,symbol,refseq.rna'):
    if gene_symbol:
     return mg.query('symbol:%s' % gene_symbol, species='human')
    elif gene_id:
     return mg.getgene(gene_id, values)

def vep_anno(chrom, pos, ref, alt):
    # doing vep annotation
    # costructing vep keys
    # and a v_id => vep_key dict
    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    print 'vep annotation'
    vep_keys_dict = {}
    vep_keys_list = []
    vep_key = '{} {} . {} {} . . .'.format(chrom, pos, ref, alt)
    #vep_keys_dict[variant] = vep_key
    #vep_keys_list.append(vep_key)
    # construct query
    query = '{ "variants" : ["' + '", "'.join(vep_key) + '" ] }'
    # do the query
    r = requests.post(server+ext, headers=headers, data=query)
    # check the query
    if not r.ok:
        r.raise_for_status()
        return
    # parse VEP
    # has potential more than one transcripts.
    #  pick the first and most damaging one
    #first_transcript = g['VEP']['transcript_consequences'][0]
    #this['consequence'] = [g['VEP']['most_severe_consequence'], {'impact': first_transcript.get('impact', 'NA'), 'polyphen': first_transcript.get('polyphen_prediction','NA'), 'sift': first_transcript.get('sift_prediction','NA')}]
    #aa = first_transcript.get('amino_acids','')
    #if aa:
       #this['p_change'] = {'pvar': 'p.%s %s' % (first_transcript['protein_start'], aa), 'cvar': 'c.%s %s' % (first_transcript['cds_start'], first_transcript['codons']), 'gene_id': gene_id, 'transcript_id': first_transcript['transcript_id']}
    #else:
       # not a aa change
       #this['p_change'] = {'pvar': 'NA', 'gene_id': gene_id, 'transcript_id':first_transcript['transcript_id']}
    # has default, weak vep_annotations
    #vep = g['vep_annotations'][0]
    #this['consequence'] = [vep['Consequence'], {}]
    #aa = vep['Amino_acids']
    #this['p_change'] = {'pvar': 'p.%s %s' % (vep['Protein_position'], aa), 'cvar': 'c.%s %s' % (vep['CDS_position'], vep['Codons']), 'gene_id': gene_id, 'transcript_id': vep['Feature']}
    #if debug: print(r.json())
    return r.json()
    #else: return json.loads(r.json())


def rvs_anno(chrom, pos, ref, alt):
    # annotate variants with rvs, and add the emtpy records to bulk_vep
    ##################
    # construct vkeys for rvs
    all_vkeys_list = []
    # dict for v_id => rvs_id
    rvs_dict = {}
    # get the functional arguments of v2k function
    rvs_id=VarCharKey.v2k(chrom, int(pos), int(pos), alt)
    all_vkeys_list.append(rvs_id)
    print all_vkeys_list
    all_vkeys = ','.join(all_vkeys_list)
    ##################
    # request rvs
    rvs_vars = {}
    for mode in ['impact','prediction']:
        url='https://rvs.u.hpc.mssm.edu/rest/{}/vkey/{}'.format(mode,all_vkeys)
        print url
        r = requests.get(url, headers={ "Content-Type" : "application/json"})
        rvs_vars[mode]=r.json()
    # parse RVS record
    #rvs=rvs[0]
    #VAR['p_change']={'gene_id':rvs['gene_id'],'cvar':rvs['hgvs_c'],'pvar':rvs['hgvs_p'],'transcript_id':rvs['enst']}
    #VAR['consequence']=rvs['effect']
    #print 'NO RVS, USING CSV INFO'
    #ENSG00000164256:ENST00000296682:exon11:c.2497_2580del:p.833_860del
    #variant=dict()
    #variant['p_change'] = {'pvar': g['RVS']['impact']['hgvs_p'], 'cvar': g['RVS']['impact']['hgvs_c'], 'gene_id': gene_id, 'transcript_id': g['RVS']['impact']['enst']}
    #variant['consequence'] = [g['RVS']['impact']['effect'], {'impact': g['RVS']['impact']['impact']} ]
    #if g['RVS'].get('prediction',{}): temp = g['RVS']['prediction'] 
    #for pred in ['Polyphen2_HDIV', 'SIFT', 'CADD', 'MutationTaster', 'ensemble_prediction', 'FATHMM', 'MutationAssessor', 'phastCons', 'GWAVA_region', 'Polyphen2_HVAR']: this['consequence'][1][pred] = temp[pred]
    return rvs_vars


def exac_anno(var,update=True):
    #if v and 'EXAC' in v and 'allele_freq' in v['EXAC']: return v['EXAC']
    attempt = 5
    while attempt:
        try:
            r=requests.get('http://exac.broadinstitute.org/variant/%s'%var)
            break
        except requests.ConnectionError:
            print 'query exac connectionError, retry'
            attempt -= 1
            time.sleep(1)
    if not r: raise 'Too many attempts to query exac, fail at exac_anno'
    #http://exac.hms.harvard.edu/rest/variant/1-1271580-C-T
    p=r.text
    VAR={}
    m=re.compile('This variant is not found in ExAC.').search(p)
    if m:
        print var, 'not in exac'
        return VAR
    else:
        print var, 'in exac'
    m=re.compile("window.variant\s*=\s*({.*})\s*;").search(p)
    if m: variant=json.loads(m.group(1))
    #VAR['variant']=variant
    VAR['genes']=variant['genes']
    VAR['vep_annotations']=variant['vep_annotations']
    VAR['pop_homs']=variant['pop_homs']
    VAR['pop_acs']=variant['pop_acs']
    VAR['pop_ans']=variant['pop_ans']
    VAR['pop_af']={}
    for k in VAR['pop_homs']:
        try:
            VAR['pop_af'][k]=float(VAR['pop_acs'][k])/float(VAR['pop_ans'][k])
        except:
            VAR['pop_af'][k]=None
    VAR['total_homs']=sum([int(VAR['pop_homs'][k]) for k in VAR['pop_homs']])
    VAR['total_acs']=sum([int(VAR['pop_acs'][k]) for k in VAR['pop_acs']])
    VAR['total_ans']=sum([int(VAR['pop_ans'][k]) for k in VAR['pop_ans']])
    if float(VAR['total_ans'])==0: 
        VAR['allele_freq']=None
    else:
        VAR['allele_freq']=float(VAR['total_acs'])/float(VAR['total_ans'])
    m=re.compile("window.consequence\s*=\s*({.*})\s*;").search(p)
    if m: csq=json.loads(m.group(1))
    #VAR['consequence']=csq
    m=re.compile("window.metrics\s*=\s*({.*})\s*;").search(p)
    if m: metrics=json.loads(m.group(1))
    #VAR['metrics']=metrics
    return VAR


