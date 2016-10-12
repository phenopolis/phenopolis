
csq_order = ["transcript_ablation",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "initiator_codon_variant",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "non_coding_exon_variant",  # deprecated
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "nc_transcript_variant",  # deprecated
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "regulatory_region_variant",
        "feature_elongation",
        "feature_truncation",
        "intergenic_variant",
        "start_lost",
        'protein_altering_variant']

csq_order=[
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'start_lost',
    'protein_altering_variant' ]

from collections import Counter
import pymongo

conn = pymongo.MongoClient(host='localhost', port=27017)
db=conn['uclex-old']

for csq in csq_order:
    print(csq, db.variants.find({'vep_annotations.Consequence':csq,'filter':'PASS'}).count(),)

for csq in csq_order:
    print(csq, db.variants.find({'vep_annotations.Consequence':csq,'in_exac':True,'filter':'PASS'}).count(),)

for csq in csq_order:
    print(csq, db.variants.find({'vep_annotations.Consequence':csq,'in_exac':False,'filter':'PASS'}).count(),)



ptvs=[
'stop_gained',
'frameshift_variant',
'stop_lost',
'inframe_insertion',
'inframe_deletion',
#'missense_variant',
'splice_region_variant',
'incomplete_terminal_codon_variant',
'stop_retained_variant',
#'synonymous_variant',
#'coding_sequence_variant',
#'start_lost',
#'protein_altering_variant',
]



print(db.variants.find({'vep_annotations.Consequence': {'$in': ptvs},'filter':'PASS','in_exac':False}).count(),)


print(db.variants.find({'vep_annotations.Consequence': 'stop_gained'}).count(),)

ptv_genes=[v['genes'][0] for v in db.variants.find( {'vep_annotations.Consequence': 'frameshift_variant','filter':'PASS','in_exac':False})]

gene_counts=Counter( ptv_genes )

import operator
gene_counts_sorted=sorted(gene_counts.items(),key=operator.itemgetter(1),reverse=True)

for g, c, in gene_counts_sorted:
    print(db.genes.find_one({'gene_id':g})['gene_name_upper'],c,)



#db.variants.aggregate( [ {'$match':{'vep_annotations.Consequence': {'$in': [ 'stop_gained', 'frameshift_variant', 'stop_lost', 'inframe_insertion', 'inframe_deletion', 'splice_region_variant', 'incomplete_terminal_codon_variant', 'stop_retained_variant' ]},'filter':'PASS'}}, {'$out':'protein_truncating_variants'}] )

import pysam
for v in db.non_exac_ptv.find():
    # who has the variant?
    g=v['genes'][0]
    gene=db.genes.find_one({'gene_id':g})['gene_name_upper']
    variant_str=v['variant_id']
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    g=[r.split('\t') for r in records]
    if len(g)==0: continue
    geno=dict(zip(headers, g[0]))
    WT=[h for h in geno if geno[h].split(':')[0]=='0/0']
    HETS=[h for h in geno if geno[h].split(':')[0]=='0/1']
    HOMS=[h for h in geno if geno[h].split(':')[0]=='1/1']
    if (len(HOMS)>len(WT)): db.variants.remove({'variant_id':variant_str})
    features=dict()
    for s in samples:
        f=db.patients.find_one({'external_id':s})
        if f is None:
            continue
        if 'features' not in f: features[s]=['HP:000001']
        else: features[s]=[x['id'] for x in f['features']]
    Counter()
    # get the hpo dissimilarity between all individuals
    print(variant_str,gene,len(WT),len(HETS),len(HOMS))



individual_counts={}


import pysam
for v in db.non_exac_ptv.find():
    # who has the variant?
    gene_id=v['genes'][0]
    gene=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
    variant_str=v['variant_id']
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    g=[r.split('\t') for r in records]
    if len(g)==0: continue
    geno=dict(zip(headers, g[0]))
    WT=[h for h in geno if geno[h].split(':')[0]=='0/0']
    HETS=[h for h in geno if geno[h].split(':')[0]=='0/1']
    HOMS=[h for h in geno if geno[h].split(':')[0]=='1/1']
    if (len(HOMS)>len(WT)): continue
    for h in HETS+HOMS:
        individual_counts[h]=individual_counts.get(h,0)+1


individual_gene_counts={}
import pysam
for v in db.non_exac_ptv.find():
    # who has the variant?
    g=v['genes'][0]
    gene=db.genes.find_one({'gene_id':g})['gene_name_upper']
    variant_str=v['variant_id']
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    g=[r.split('\t') for r in records]
    if len(g)==0: continue
    geno=dict(zip(headers, g[0]))
    WT=[h for h in geno if geno[h].split(':')[0]=='0/0']
    HETS=[h for h in geno if geno[h].split(':')[0]=='0/1']
    HOMS=[h for h in geno if geno[h].split(':')[0]=='1/1']
    if (len(HOMS)>len(WT)): continue
    for h in HETS+HOMS:
        individual_gene_counts[h]=individual_gene_counts.get(h,[])+[gene]

individual_gene_counts2=individual_gene_counts

individual_gene_counts3=dict([(i,Counter(individual_gene_counts[i])) for i in individual_gene_counts])

import csv
writer = csv.DictWriter(file('nonexac_ptv_genes.csv','wb'), fieldnames=['individual','gene','count'])
writer.writeheader()
#outfile  = csv.writer(file('nonexac_ptv_genes.csv','wb'), delimiter=',')
for i in individual_gene_counts3:
    for j in (individual_gene_counts3[i]):
        writer.writerow({'individual':i,'gene':j,'count':individual_gene_counts3[i][j]})


exac_genes=[g['genes'][0] for g in db.protein_truncating_variants.find({in_exac:True},{'genes':1,_id:False})]
nonexac_genes=[g['genes'][0] for g in db.protein_truncating_variants.find({in_exac:False},{'genes':1,_id:False})]





