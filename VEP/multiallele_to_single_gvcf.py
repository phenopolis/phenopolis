### Hopefully vcfbreakmulti, from https://github.com/ekg/vcflib, can take care of this now.
#! /bin/env python
from __future__ import print_function
import sys
import argparse

#these 9 column headers are standard to all VCF files
#the samples headers depend on the number of samples in the file
#which we find out once we read the #CHROM line
SAMPLE_HEADERS=[]

parser=argparse.ArgumentParser(description='Arguments to multiallele_to_single_gvcf.py')
parser.add_argument('--GQ',default=None,type=int)
parser.add_argument('--DP',default=None,type=int)
parser.add_argument('--headers',default=','.join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']),type=str)

args=parser.parse_args()
#args.GQ
#args.DP
#STD_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
STD_HEADERS=args.headers.strip().split(',')

def print_line(s):
    #if s['ID']=='.':
    s['ID']='-'.join([s['CHROM'],s['POS'],s['REF'][:4],s['ALT'][:4]])
    print( *([s[h] for h in STD_HEADERS+SAMPLE_HEADERS]), sep='\t' )

print('##fileformat=VCFv4.1')
for line in sys.stdin:
    line=line.strip()
    #remove beginning
    #lines which start with '###'
    #are not tab separated
    if line.startswith('##'):
        #print(line)
        continue
    #this is tab separated line
    s=line.split("\t")
    #header line yay!: #CHROM ...
    if line.startswith('#'):
        print(line)
        headers=s
        headers[0]=headers[0].strip('#')
        #the first 9 names in the header are standard (see above)
        if (headers[0:len(STD_HEADERS)] != STD_HEADERS):
            print( headers[0:len(STD_HEADERS)],  STD_HEADERS )
            raise 'hell'
        #everything else in the header is a sample name
        SAMPLE_HEADERS=headers[len(STD_HEADERS):]
        continue
    #you can now access elements by their header name
    s=dict(zip(headers,s))
    #split alternate alleles
    alternative_alleles=s['ALT'].split(',')
    #I would expect each sample to be formatted according to the format column.
    #However I found in practise this not true.
    #Only the first two fields: GT (genotype) and AD (allele depth)
    #are compulsory, the remaining fields, including DP (total depth),
    #can be missing.
    #So I calculate DP from AD instead.
    #I ignore the other fields for now.
    #Actually, if the genotype is missing then none of the fields are compulsory.
    #In fact AD can be '.' even when geno is 0/0 if using HaplotypeCaller with multiple BAMs.
    for h in SAMPLE_HEADERS:
        d=dict(zip(s['FORMAT'].split(':'),s[h].split(':')))
        GT=d['GT'].replace('|','/')
        AD=d.get('AD',','.join(['0']*(len(alternative_alleles)+1)))
        if 'DP' in d:
            DP=(d['DP'])
        else:
            DP=str(sum([int(x) if x!='.' else 0 for x in AD.split(',')]))
        GQ=d.get('GQ','')
        # if GQ is '.' the must be missing
        if GQ == '.' and GT != './.': raise 'hell'
        # if fails QC set to missing
        if ((args.DP and int(DP) < args.DP) or (args.GQ and GQ!='.' and int(GQ) < args.GQ)): GT='./.'
        PL=d.get('PL','')
        s[h]=':'.join([GT,AD,DP,GQ,PL])
    #if single alternate allele
    #then just print out as normal
    if len(alternative_alleles)<=1:
        s['FORMAT']='GT:AD:DP:GQ:PL'
        print_line(s)
        continue
    #otherwise split over as many lines as there are alternative alleles
    #each genotype then takes a 1 if matches the alternative or a 0 otherwise
    alleles=[s['REF']]+alternative_alleles
    n2geno=dict(zip(['.']+[str(i) for i in range(0,len(alleles))],['.']+alleles))
    for idx, alt in enumerate(alternative_alleles):
        s['ALT']=alt
        #recode GT
        s1=s.copy()
        for h in SAMPLE_HEADERS:
            d=dict(zip(s1['FORMAT'].split(':'),s1[h].split(":")))
            #length of allele depth is 2 where first is always REF allele depth
            #and second can be either ALT
            if 'AD' in d and d['AD']!='' and d['AD']!='.':
                AD=d['AD'].split(',')
                AD[1]=AD[idx+1]
                AD=','.join(AD[:2])
            else:
                AD='0'
            # matches REF -> 0
            # matches ALT -> 1
            # matches anything else -> .
            #sorted so that "." precedes number
            GT='/'.join(sorted([{s['REF']:'0',s['ALT']:'1'}.get(n2geno[g],'.') for g in d['GT'].split('/')]))
            if AD!='.':
                DP=str(sum([int(x) for x in AD.split(',')]))
            else:
                DP=0
            s1[h]=':'.join( [GT, AD, DP] )
            s1['FORMAT']='GT:AD:DP'
        print_line(s1)





