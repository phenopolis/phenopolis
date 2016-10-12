
import pysam

# VCF query
def vcf_query(chrom=None, pos=None, ref=None, alt=None, variant_str=None, individual=None, verbose=False, limit=100, release='mainset_July2016'):
    if variant_str:
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
    #tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    tb=pysam.TabixFile('/slms/gee/research/vyplab/UCLex/%s/%s_chr%s.vcf.gz' % (release, release, chrom,))
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
            print 'POS', POS
            print 'REF', REF
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            if verbose: print 'ALT', ALT
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



# VCF query
def vcf_query2(chrom=None, pos=None, ref=None, alt=None, variant_str=None, individual=None, verbose=False, limit=100):
    if variant_str:
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
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
            print 'POS', POS
            print 'REF', REF
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            if verbose: print 'ALT', ALT
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


# VCF query
def vcf_query3(chrom=None, pos=None, ref=None, alt=None, variant_str=None, individual=None, verbose=False, limit=100):
    if variant_str:
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
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
            print 'POS', POS
            print 'REF', REF
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            if verbose: print 'ALT', ALT
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


def vcf_query_gene():
    tb=pysam.TabixFile('/slms/gee/research/vyplab/UCLex/%s/%s_chr%s.vcf.gz' % (RELEASE, RELEASE, gene.chrom,))
    region ='%s:%s-%s' % (str(gene.chrom), str(gene.start), str(gene.stop),)
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip('#').strip().split('\t')
    records=[dict(zip(headers,r.strip().split('\t'))) for r in tb.fetch(region)]
    print(len(records))
    records=dict([('%s-%s-%s-%s' % (r['CHROM'], r['POS'], r['REF'], r['ALT'],),r,) for r in records])

