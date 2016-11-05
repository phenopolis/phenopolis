
import sys
import re
import traceback
#from utils import *
import copy
import sys
from collections import Counter
from collections import OrderedDict
import pysam
import itertools
import psycopg2

CHROMOSOMES = ['chr%s' % x for x in range(1, 23)]
CHROMOSOMES.extend(['chrX', 'chrY', 'chrM'])
CHROMOSOME_TO_CODE = { item: i+1 for i, item in enumerate(CHROMOSOMES) }

def get_single_location(chrom, pos):
    """
    Gets a single location from chromosome and position
    chr must be actual chromosme code (chrY) and pos must be integer

    Borrowed from xbrowse
    """
    return CHROMOSOME_TO_CODE[chrom] * int(1e9) + pos

def get_xpos(chrom, pos):
    """
    Borrowed from xbrowse
    """
    if not chrom.startswith('chr'):
        chrom = 'chr{}'.format(chrom)
    return get_single_location(chrom, int(pos))

def get_minimal_representation(pos, ref, alt): 
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets, 
    with different combinations of alternate alleles, can always be matched directly. 

    Note that chromosome is ignored here - in xbrowse, we'll probably be dealing with 1D coordinates 
    Args: 
        pos (int): genomic position in a chromosome (1-based)
        ref (str): ref allele string
        alt (str): alt allele string
    Returns: 
        tuple: (pos, ref, alt) of remapped coordinate
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1: 
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt 
#cur.execute("CREATE TABLE variants_patients (id serial PRIMARY KEY, variant_id integer, data varchar);")
#cur.execute("SELECT * FROM test;") 

#chrom=sys.argv[1]


conn = psycopg2.connect(database="uclex")
cur = conn.cursor()

chrom=sys.argv[1]
print(chrom)

# uclex files
vcf_reader = pysam.VariantFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom)
for v in vcf_reader:
    for alt in v.alts:
        variant=dict()
        variant['variant_id']='%s-%s-%s-%s' % (v.chrom, v.pos,v.ref,alt,)
        GT=[v.samples[s]['GT'] for s in v.samples]
        miss=[s for s in v.samples if v.samples[s]['GT'].count(None)==2]
        variant['miss_count']=len(miss)
        wt=[s for s in v.samples if v.samples[s]['GT'].count(0)==2]
        variant['wt_count']=len(wt)
        het=[s for s in v.samples if v.samples[s]['GT'].count(1)==1]
        variant['het_count']=len(het)
        hom=[s for s in v.samples if v.samples[s]['GT'].count(1)==2]
        variant['hom_count']=len(hom)
        #variant['MUT']=variant['HET']+variant['HOM']
        # Make a copy of the info_field dict - so all the original data remains
        # Add some new keys that are allele-specific
        #pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)
        variant['chrom'] = v.chrom
        variant['pos'] = v.pos
        variant['rsid'] = v.rid
        variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
        variant['ref'] = v.ref
        variant['alt'] = alt
        variant['xstart'] = variant['xpos']
        variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
        variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
        variant['orig_alt_alleles'] = ','.join([ '{}-{}-{}-{}'.format(variant['chrom'], *get_minimal_representation(v.pos, v.ref, x)) for x in v.alts ])
        variant['filter'] = ','.join([k for k in v.filter.keys() if k.strip()])
        variant['site_quality'] = float(v.qual)
        print(variant['variant_id'])
        cur.execute( """ INSERT INTO variants (variant_id, chrom, pos, xstart, xstop, ref, alt, filter, site_quality, miss_count, wt_count, het_count, hom_count) VALUES (%(variant_id)s, %(chrom)s, %(pos)s, %(xstart)s, %(xstop)s, %(ref)s, %(alt)s, %(filter)s, %(site_quality)s, %(miss_count)s, %(wt_count)s, %(het_count)s, %(hom_count)s) """, variant)
        conn.commit()
        if variant['wt_count']<variant['hom_count']:
            print('flipped variant')
            continue
        for p in het:
            cur.execute(""" INSERT INTO het_variants_patients (variant_id, patient_id) VALUES (%s, %s) """, (variant['variant_id'], p))
        for p in hom:
            cur.execute(""" INSERT INTO hom_variants_patients (variant_id, patient_id) VALUES (%s, %s) """, (variant['variant_id'], p))
        conn.commit()

cur.close()
conn.close()

