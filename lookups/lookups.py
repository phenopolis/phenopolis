import re
#from utils import *
import itertools
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import csv
#hpo lookup
import phizz
import random
import pickle
import hashlib
import pprint
import orm

SEARCH_LIMIT = 10000
# massive genes?
#UNSUPPORTED_QUERIES = ['TTN', 'ENSG00000155657', 'CMD1G', 'CMH9', 'CMPD4', 'FLJ32040', 'LGMD2J', 'MYLK5', 'TMD', u'ENST00000342175', u'ENST00000359218', u'ENST00000342992', u'ENST00000460472', u'ENST00000589042', u'ENST00000591111']

def lookup_patient(db,user,external_id):
    external_ids=db.users.find_one({'user':user},{'external_ids':1})['external_ids']
    return external_id in external_ids

def xpos_to_pos(xpos): return int(xpos % 1e9)

def get_gene(db, gene_id):
    print(gene_id)
    for g in db.genes.find({'gene_id': gene_id}): print(g)
    #return g
    return db.genes.find_one({'gene_id': gene_id}, projection={'_id': False})


def get_gene_by_name(db, gene_name):
    # try gene_name field first
    gene = db.genes.find_one({'gene_name': gene_name}, projection={'_id': False})
    if gene: return gene
    # if not, try gene['other_names']
    return db.genes.find_one({'other_names': gene_name}, projection={'_id': False})


def get_transcript(db, transcript_id):
    transcript = db.transcripts.find_one({'transcript_id': transcript_id}, projection={'_id': False})
    if not transcript:
        return None
    transcript['exons'] = get_exons_in_transcript(db, transcript_id)
    return transcript


def get_raw_variant(db, xpos, ref, alt, get_id=False):
    return db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection={'_id': get_id})

def get_variant(db, variant_id):
    return db.variants.find_one({'variant_id':variant_id})


def get_variant(db, xpos, ref, alt):
    variant = get_raw_variant(db, xpos, ref, alt, False)
    print(variant)
    if variant is None or 'rsid' not in variant: return variant
    if variant['rsid'] == '.' or variant['rsid'] is None:
        rsid = db.dbsnp.find_one({'xpos': xpos})
        if rsid:
            variant['rsid'] = 'rs%s' % rsid['rsid']
    return variant

def get_variants_from_dbsnp(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        rsid = int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    position = db.dbsnp.find_one({'rsid': rsid})
    if position:
        variants = list(db.variants.find({'xpos': {'$lte': position['xpos'], '$gte': position['xpos']}}, projection={'_id': False}))
        if variants:
            #add_consequence_to_variants(variants)
            return variants
    return []


def get_coverage_for_bases(db, xstart, xstop=None):
    """
    Get the coverage for the list of bases given by xstart->xstop, inclusive
    Returns list of coverage dicts
    xstop can be None if just one base, but you'll still get back a list
    """
    if xstop is None:
        xstop = xstart
    coverages = {
        doc['xpos']: doc for doc in db.base_coverage.find(
            {'xpos': {'$gte': xstart, '$lte': xstop}},
            projection={'_id': False}
        )
    }
    ret = []
    for i in range(xstart, xstop+1):
        if i in coverages:
            ret.append(coverages[i])
        else:
            ret.append({'xpos': i, 'pos': xpos_to_pos(i)})
    for item in ret:
        item['has_coverage'] = 'mean' in item
        del item['xpos']
    print '+++++++++++++++++++++++++++'
    temp = db.base_coverage.find({'xpos': {'$gte': xstart, '$lte': xstop}})
    from bson.json_util import dumps
    dumps(temp)
    print xstart
    print xstop
    print '+++++++++++++++++++++++++++++'
    return ret


def get_coverage_for_transcript(db, xstart, xstop=None):
    """
    :param db:
    :param genomic_coord_to_exon:
    :param xstart:
    :param xstop:
    :return:
    """
    coverage_array = get_coverage_for_bases(db, xstart, xstop)
    # only return coverages that have coverage (if that makes any sense?)
    # return coverage_array
    #print '+++++++++++++++++++++++++'
    #print coverage_array
    #print '+++++++++++++++++++++++++'
    covered = [c for c in coverage_array if c['has_coverage']]
    for c in covered: del c['has_coverage']
    return covered


def get_constraint_for_transcript(db, transcript):
    return db.constraint.find_one({'transcript': transcript}, projection={'_id': False})


def get_awesomebar_suggestions(g, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = (r for r in g.autocomplete_strings if regex.match(r))
    results = itertools.islice(results, 0, 20)
    #db.hpo.find({'name':/.*retinal dystrophy.*/})
    return list(results)


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
R2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
R3 = re.compile(r'^(\d+|X|Y|M|MT)$')
R4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')



def get_genes_in_region(db, chrom, start, stop):
    """
    Genes that overlap a region
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    genes = db.genes.find({ 'xstart': {'$lte': xstop}, 'xstop': {'$gte': xstart}, }, projection={'_id': False})
    return list(genes)


def get_variants_in_region(db, chrom, start, stop):
    """
    Variants that overlap a region
    Unclear if this will include CNVs
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    variants = list(db.variants.find({ 'xpos': {'$lte': xstop, '$gte': xstart}
    }, projection={'_id': False}, limit=SEARCH_LIMIT))
    #add_consequence_to_variants(variants)
    return list(variants)



def remove_extraneous_information(variant):
    return
    del variant['genotype_depths']
    del variant['genotype_qualities']
    del variant['transcripts']
    del variant['genes']
    del variant['orig_alt_alleles']
    del variant['xpos']
    del variant['xstart']
    del variant['xstop']
    del variant['site_quality']
    del variant['vep_annotations']



def get_transcripts_in_gene(db, gene_id):
    """
    """
    return list(db.transcripts.find({'gene_id': gene_id}, projection={'_id': False}))


def get_exons_in_transcript(db, transcript_id):
    # return sorted(
    #     [x for x in
    #      db.exons.find({'transcript_id': transcript_id}, fields={'_id': False})
    #      if x['feature_type'] != 'exon'],
    #     key=lambda k: k['start'])
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, fields={'_id': False})), key=lambda k: k['start'])


def get_hpo_patients(hpo_db, patients_db, hpo_id, cached=True,verbose=False):
    """
    Get patients with HPO term.
    """
    if cached:
        return [p for p in patients_db.patients.find({'external_id':{'$in':patients_db.hpo_cache.find_one({'hpo_id':hpo_id})['external_id']}}) if 'external_id' in p]
    if 'HP:0000001' == hpo_id: return [p for p in patients_db.patients.find() if 'external_id' in p]
    patients = [p for p in patients_db.patients.find({'features.id':hpo_id}) for f in p['features']  if f['id']== hpo_id and f['observed']=='yes']
    if verbose: print(hpo_id,len(patients))
    for r in hpo_db.hpo.find({'is_a':hpo_id}):
        for i in r['id']: patients+=list(itertools.chain(get_hpo_patients(hpo_db,patients_db,i,cached=cached,verbose=verbose))) 
    #remove duplicates
    patients={v['external_id']:v for v in patients}.values()
    return patients

# return hpo terms found in people in which variant is found
def get_hpo(variant_str):
    samples=get_samples(variant_str)
    #chrom,pos,ref,alt,=str(variant_str.strip()).split('-')
    d=csv.DictReader(file('/data/uclex_data/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
    hpo=[]
    for r in d:
        if r['sample'] not in samples: continue
        pheno=r['phenotype']
        print((r['sample'],pheno,))
        if pheno.startswith('HP'):
            hpo+=[phizz.query_hpo([pheno])]
        elif pheno.startswith('MIM'):
            hpo+=[phizz.query_disease([pheno])]
    return(hpo)

def get_hpo_children(hpo_db, hpo_id):
    hpo=[hpo_db.hpo.find_one({'id':hpo_id})]
    for r in hpo_db.hpo.find({'is_a':hpo_id}):
        for i in r['id']:
            hpo+=list(itertools.chain(get_hpo_children(hpo_db,i))) 
    #remove duplicates
    hpo={h['id'][0]:h for h in hpo}.values()
    return hpo

def replace_hpo(hpo_db, hpo):
    # some hpo_ids are obsolete.
    record = hpo_db.hpo.find_one({'id':hpo[0]})
    if not record:
        print 'no record in replace_hpo'
        print hpo
    if 'replaced_by' in record:
        new = hpo_db.hpo.find_one({'id':record['replaced_by'][0]})
        return [new['id'][0], new['name'][0]]
    else:
        return hpo

def get_hpo_ancestors(hpo_db, hpo_id):
    """
    Get HPO terms higher up in the hierarchy.
    """
    h=hpo_db.hpo.find_one({'id':hpo_id})
    #print(hpo_id,h)
    if 'replaced_by' in h:
        # not primary id, replace with primary id and try again
        h = hpo_db.hpo.find_one({'id':h['replaced_by'][0]})
    hpo=[h]
    if 'is_a' not in h: return hpo
    for hpo_parent_id in h['is_a']:
        #p=hpo_db.hpo.find({'id':hpo_parent_id}):
        hpo+=list(itertools.chain(get_hpo_ancestors(hpo_db,hpo_parent_id))) 
    #remove duplicates
    hpo={h['id'][0]:h for h in hpo}.values()
    return hpo

def get_hpo_ancestors_array(hpo_db, hpo_id):
    # return an array of ids, instead of array of dicts
    anc = get_hpo_ancestors(hpo_db, hpo_id)
    result = []
    for a in anc:
        result.extend(a['id'])
    return result

def get_hpo_size_freq(freq_file):
    # read freq file
    # result = {'HP:0000345':{size: 456, freq: 0.1, raw: 456/4500}}
    hpo_freq = {}
    inf = open(freq_file, 'r')
    for l in inf:
        l = l.rstrip().split('\t')
        nums = l[1].split('/')
        size = int(nums[0])
        tot = float(nums[1])
        hpo_freq[l[0]] = {'size': size, 'freq': size/tot, 'raw': l[1]}
    return hpo_freq

def get_hpo_common_ancestors(hpo_db, h1, h2):
    # return a list of hpo ids for h1 and h2's common ancestors
    a1 = get_hpo_ancestors(hpo_db, h1)
    a2 = get_hpo_ancestors(hpo_db,h2)
    an1 = []
    an2 = []
    for a in a1:
        an1.extend(a['id'])
    for a in a2:
        an2.extend(a['id'])
    return list(set(an1) & set(an2))

def get_hpo_nearest_common_ancestors(hpo_db, h1, h2, hpo_freq):
    # given hpo_freq, find out a list of nearest common ancestors
    common_ans = get_hpo_common_ancestors(hpo_db, h1, h2)
    freqs = [hpo_freq[h] for h in common_ans]
    min_freq = min(freqs)
    inds = [i for i, v in enumerate(freqs) if v == min_freq]
    return [common_ans[i] for i in inds]

def hpo_minimum_set(hpo_db, hpo_ids=[]):
    '''
    minimize the hpo sets
    results = {'HP:0000505': [ancestors]}
    '''
    hpo_ids = list(set(hpo_ids))
    results = dict([(hpo_id, [ h['id'][0] for h in get_hpo_ancestors(hpo_db, hpo_id)],) for hpo_id in hpo_ids])
    # minimise
    bad_ids = []
    for i in range(len(hpo_ids)):
        for j in range(i+1,len(hpo_ids)):
            if hpo_ids[i] in results[hpo_ids[j]]:
                # i is j's ancestor, remove
                bad_ids.append(hpo_ids[i])
                break
            if hpo_ids[j] in results[hpo_ids[i]]:
                # j is i's ancestor, remove
                bad_ids.append(hpo_ids[j])
    return list(set(hpo_ids) - set(bad_ids))


def get_patient_hpo(hpo_db,patients_db, patient_id,ancestors=True):
    """
    Get complete hierarchy of HPO terms for patient.
    """
    p=patients_db.patients.find_one({'external_id':patient_id})
    if 'features' not in p: return []
    if ancestors:
        hpo_ancestors=[]
        for hpo_ids in [f['id'] for f in p['features'] if f['observed']=='yes']:
            hpo_ancestors+=get_hpo_ancestors(hpo_db,hpo_ids)
        # remove duplicates
        hpo_ancestors={h['id'][0]:h for h in hpo_ancestors}.values()
        return hpo_ancestors
    else:
        return [ hpo_db.hpo.find_one({'id':f['id']}) for f in p['features'] if f['observed']=='yes']

def get_gene_hpo(hpo_db,gene_name,dot=True):
    """
    Get all HPO terms linked to gene name, including ancestors.
    and return as dot string for plotting if dot is True.
    """
    hpo_ids=[hpo['HPO-Term-ID'] for hpo in hpo_db.OMIM_ALL_FREQUENCIES_genes_to_phenotype.find({'entrez-gene-symbol':gene_name})]
    if not hpo_ids:
        hpo_ids=hpo_db.genes_pheno.find_one({'gene':gene_name})
        # no hpo linked to gene
        if hpo_ids is None: hpo_ids=[]
        else: hpo_ids=hpo_ids['hpo']
    hpo_ancestors=[get_hpo_ancestors(hpo_db,hid) for hid in hpo_ids]
    hpo_ancestors=list(itertools.chain(*hpo_ancestors)) 
    # remove duplicates
    hpo_ancestors={h['id'][0]:h for h in hpo_ancestors}.values()
    hpo_string="digraph {"
    for h in hpo_ancestors:
        hpo_id=h['id'][0]
        hpo_label=h['name'][0]
        #hpo_count=0
        hpo_string+= '"{}" [style="filled", fixedsize="true", fontsize="15", shape="circle", width="0.75", fillcolor="powderblue", label="{}\n{}", color="transparent"];\n'.format(hpo_id,hpo_label,hpo_id)
    for h in hpo_ancestors:
        hpo_id=h['id'][0]
        if 'is_a' not in h: continue
        for anc in h['is_a']:
            hpo_string+='"{}" -> "{}" [color="#000000", lty="solid"];\n'.format(anc,hpo_id)
    hpo_string+= '}'
    if dot:
        return hpo_string
    else:
        return hpo_ancestors


# get hpo terms shared between patients
def common_hpo(hpo_db,patients_db,patient_ids):
    terms_by_patient=[get_patient_hpo(hpo_db,patients_db,pid) for pid in patient_ids]
    # intersection of lists
    common_hpo_term_ids=frozenset.intersection(*[frozenset([y['id'][0] for y in x]) for x in terms_by_patient])
    # remove ancestors
    #get_hpo_ancestors(hpo_db, hpo_id):
    # lookup hpo terms
    common_hpo_terms=[hpo_db.hpo.find_one({'id':hpo_id}) for hpo_id in common_hpo_term_ids]
    return common_hpo_terms

# get union of hpo terms seen in  patients
def union_hpo(hpo_db,patients_db,patient_ids):
    terms_by_patient=[get_patient_hpo(hpo_db,patients_db,pid) for pid in patient_ids]
    #flatten lists
    terms_by_patient=list(itertools.chain(*terms_by_patient)) 
    # intersection of lists
    terms_by_patient={h['id'][0]:h for h in terms_by_patient}.values()
    return terms_by_patient


# VCF gene query
def variants_in_gene_vcf(gene_symbol):
    import mygene
    mg = mygene.MyGeneInfo()
    g=mg.query('symbol:%s' % gene_symbol, fields='exons', species='human')
    print g
    exons=g['hits'][0]['exons']
    for transcript in exons:
        yield (transcript, exons[transcript],)

def get_patient_observed_hpo(patient, patient_db):
    # returns [('HP:0000001', 'hell yeah')]
    this_patient = patient_db.patients.find_one({'external_id':patient}) 
    result = [(None, None)]
    if not this_patient:
        #print 'ERROR: %s not in patients db' % patient
        pass
    else:
        if 'features' not in this_patient:
            print 'WARNING: features not in ' + patient
        p_features = this_patient.get('features', [{'id':'HP:0000001', 'label':'All', 'observed': 'yes' }])
        result = [(f['id'], f['label']) for f in p_features if f['observed']=='yes']
    return result



