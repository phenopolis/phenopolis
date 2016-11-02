#! /bin/env python
from __future__ import print_function
import sys
import json
import pysam
import exac


def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

# Takes the output of VEP and reformats

headers=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']


# VCF query
def vcf_query(chrom=None, pos=None, ref=None, alt=None, variant_str=None, individual=None, verbose=False, limit=100, release='mainset_July2016'):
    if variant_str:
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
    #tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    tb=pysam.TabixFile('/SAN/vyplab/UCLex/current/%s_chr%s.vcf.gz' % (release, chrom,))
    #mainset_February2016_chrX_filtered.vcf.gz
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    records=[r.split('\t') for r in records]
    def response(POS, REF, ALT, index, geno, chrom, pos):
        alleles=[geno['REF']]+geno['ALT'].split(',')
        homozygous_genotype='/'.join([str(index),str(index)])
        heterozygous_genotype='/'.join(['0',str(index)])
        variant=dict()
        variant['POS']=POS
        variant['REF']=REF
        variant['ALT']=ALT
        variant['index']=index
        variant['variant_id']='-'.join([str(chrom),str(POS),variant['REF'],variant['ALT']])
        variant['synonym_variant_id']='{}-{}-{}-{}'.format(str(chrom), str(pos), ref, alt,)
        variant['hgvs']='chr%s:g.%s%s>%s' % (str(chrom), str(POS), REF, ALT,)
        #print [geno[h].split(':')[0].split('/') for h in geno]
        variant['hom_samples']=[h for h in geno if geno[h].split(':')[0]==homozygous_genotype][0:limit]
        variant['HOM_COUNT']=len(variant['hom_samples'])
        variant['het_samples']=[h for h in geno if geno[h].split(':')[0]==heterozygous_genotype][0:limit]
        variant['HET_COUNT']=len(variant['het_samples'])
        variant['wt_samples']=[h for h in geno if geno[h].split(':')[0]=='0/0'][1:100]
        variant['WT_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='0/0'])
        variant['MISS_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='./.'])
        variant['allele_num']= 2*(variant['HOM_COUNT'] + variant['HET_COUNT']+variant['WT_COUNT'])
        variant['allele_count']=2*variant['HOM_COUNT'] + variant['HET_COUNT']
        if individual: variant['individual']=geno[individual]
        #variant['site_quality'] = variant['QUAL']
        #variant['filter'] = variant['FILTER']
        if variant['WT_COUNT']==0:
           variant['allele_freq'] = None
        else:
           variant['allele_freq'] = float(variant['HET_COUNT']+2*variant['HOM_COUNT']) / float(2*variant['WT_COUNT'])
        samples=variant['het_samples']+variant['hom_samples']
        #variant['hpo']=[p for p in get_db('patients').patients.find({'external_id':{'$in':samples}},{'_id':0,'features':1,'external_id':1})]
        return variant
    for r in records:
        geno=dict(zip(headers, r))
        POS=geno['POS']
        REF=geno['REF']
        if verbose:
            print('ERROR:', 'POS', POS)
            print('ERROR:', 'REF', REF)
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            if verbose: print('ALT', ALT)
            # insertion
            if ref=='-' and REF+alt==ALT: return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # deletion
            # replace leftmost
            elif alt=='-' and ALT==REF.replace(ref,''): return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # replace rightmost
            elif alt=='-' and ALT==REF[::-1].replace(ref[::-1], "", 1)[::-1]: return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            #
            elif alt=='-' and ref==REF and ALT=='*': return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            elif alt=='0' and ALT=='*' and ref==REF: return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            elif alt==ALT and ref==REF: return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            continue
    return []


def clean(d,field):
    for cons in d['transcript_consequences']:
        if field not in cons: continue
        for x in cons[field].split(','):
            if len(x.split(':'))!=2:
                sys.stderr.write(x)
                sys.stderr.write('\n')
                continue
            id,num,=x.split(':') 
            if len(id.split('>'))!=2:
                sys.stderr.write(id)
                sys.stderr.write('\n')
                continue
            ref,alt,=id.split('>')
            if alt!=d['ALT']: continue
            try:
                cons[field]=float(num)
            except:
                #eprint(num)
                print('ERROR:', 'not a number', num, id)

def freq_cleanup(d):
    for cons in d['transcript_consequences']:
        for k in [k for k in cons.keys() if k.startswith('exac') or k=='kaviar']:
            d[k]=cons[k]
    for cons in d['transcript_consequences']:
        for k in [k for k in cons.keys() if k.startswith('exac') or k=='kaviar']:
            del cons[k]

def go_cleanup(d):
    for cons in d['transcript_consequences']:
        if 'go' not in cons: continue
        cons['go']=cons['go'].split(',')

def canonical(d):
    for cons in d['transcript_consequences']:
        if 'canonical' not in cons: continue
        d['canonical_cadd']=d.get('canonical_cadd',[])+[cons.get('cadd','')]
        d['canonical_hgvsc']=d.get('canonical_hgvsc',[])+[cons.get('hgvsc','')]
        d['canonical_hgvsp']=d.get('canonical_hgvsp',[])+[cons.get('hgvsp','')]
        #grab the transript
        d['canonical_transcript']=d.get('canonical_transcript',[])+[cons['transcript_id']]
        d['canonical_gene_name_upper']=d.get('gene_name_upper',[])+[cons['gene_symbol'].upper()]
        d['canonical_gene_name_upper']=d['canonical_gene_name_upper'][0]

for l in sys.stdin:
    d=json.loads(l.strip())
    d.update(dict(zip(headers,d['input'].split('\t'))))
    d.update(dict([tuple(x.split('=')) for x in d['INFO'].split(';') if len(x.split('='))==2]))
    del d['INFO']
    del d['input']
    d['variant_id']='-'.join([d['CHROM'],d['POS'],d['REF'],d['ALT']])
    if ',' in d['ALT']:
        #eprint(d['variant_id']+' MULTIALLELIC')
        print('ERROR:', d['variant_id']+' MULTIALLELIC')
        continue
    if 'transcript_consequences' not in d:
        #eprint(d['variant_id']+' NOT CODING')
        print('ERROR:',d['variant_id']+' NOT CODING')
        continue
    clean(d,'cadd')
    clean(d,'kaviar')
    clean(d,'exac_nfe')
    clean(d,'exac_sas')
    clean(d,'exac_fin')
    clean(d,'exac_eas')
    clean(d,'exac_amr')
    clean(d,'exac_afr')
    clean(d,'exac_oth')
    clean(d,'exac_adj')
    clean(d,'1kg_eur')
    clean(d,'1kg_asn')
    clean(d,'1kg_amr')
    clean(d,'1kg_afr')
    d['genes']=list(set([cons['gene_id'] for cons in d['transcript_consequences']]))
    freq_cleanup(d)
    go_cleanup(d)
    canonical(d)
    d.update(vcf_query(variant_str=d['variant_id']))
    d['EXAC']=exac.exac_query(variant_str=d['variant_id'])
    # try convert str which have a decimal point to number
    for k in d:
        try:
            d[k] = float(d[k])
        except:
            continue
    print('JSON:', json.dumps(d))





