import rest
import vcf
import json
from operator import itemgetter
import pprint
import requests


# Note that this is the current as of v77 with 2 included for backwards compatibility (VEP <= 75)
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
'protein_altering_variant',
""]
csq_order_dict = dict(zip(csq_order, range(len(csq_order))))
rev_csq_order_dict = dict(zip(range(len(csq_order)), csq_order))

def compare_two_consequences(csq1, csq2):
    if csq_order_dict[worst_csq_from_csq(csq1)] < csq_order_dict[worst_csq_from_csq(csq2)]:
        return -1
    elif csq_order_dict[worst_csq_from_csq(csq1)] == csq_order_dict[worst_csq_from_csq(csq2)]:
        return 0
    return 1

def get_protein_hgvs(csq):
    """
    Takes consequence dictionary, returns proper variant formatting for synonymous variants
    """
    if '%3D' in csq['HGVSp']:
        try:
            amino_acids = ''.join([protein_letters_1to3[x] for x in csq['Amino_acids']])
            return "p." + amino_acids + csq['Protein_position'] + amino_acids
        except Exception, e:
            print 'Could not create HGVS for: %s' % csq
    return csq['HGVSp'].split(':')[-1]


def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst annotation (In this case, index of 'frameshift_variant', so 4)
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_index(csqs.split('&'))
    :param annnotation:
    :return most_severe_consequence_index:
    """
    return min([csq_order_dict[ann] for ann in csq_list])


def worst_csq_from_list(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return the worst annotation (In this case, 'frameshift_variant')
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_from_list(csqs.split('&'))
    :param annnotation:
    :return most_severe_consequence:
    """
    return rev_csq_order_dict[worst_csq_index(csq_list)]

def worst_csq_from_csq(csq):
    """
    Input possibly &-filled csq string (e.g. 'non_coding_exon_variant&nc_transcript_variant')
    Return the worst annotation (In this case, 'non_coding_exon_variant')
    :param consequence:
    :return most_severe_consequence:
    """
    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]


def order_vep_by_csq(annotation_list):
    print('ANNOTATION LIST',annotation_list)
    output = sorted(annotation_list, cmp=lambda x, y: compare_two_consequences(x, y), key=itemgetter('consequence_terms'))
    for ann in output:
        ann['major_consequence'] = worst_csq_from_csq(ann['consequence_terms'])
    return output



def compare_two_consequences(csq1, csq2):
    if csq_order_dict[worst_csq_from_csq(csq1)] < csq_order_dict[worst_csq_from_csq(csq2)]:
        return -1
    elif csq_order_dict[worst_csq_from_csq(csq1)] == csq_order_dict[worst_csq_from_csq(csq2)]:
        return 0
    return 1

def get_variants_by_rsid(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    variants = list([Variant(data=v) for v in db.variants.find({'rsid': rsid}, projection={'_id': False})])
    #add_consequence_to_variants(variants)
    return variants


class Variant(object):
    def __init__(self, variant_id=None, db=None,data=None):
        if variant_id is None: variant_id=data['variant_id']
        self.variant_id=str(variant_id).strip().replace('_','-')
        self.chrom, self.pos, self.ref, self.alt = variant_id.split('-')
        #q=vcf.vcf_query(variant_str=self.variant_id,)
        #if q is None: raise Exception('NOT IN VCF',self.variant_id)
        #self.__dict__.update(q)
        Variant.db=db
        data=Variant.db.variants.find_one({'variant_id':self.variant_id},projection={'_id':False})
        self.__dict__.update(data)
    def __getattribute__(self, key):
        "Emulate type_getattro() in Objects/typeobject.c"
        v = object.__getattribute__(self, key)
        if hasattr(v, '__get__'): return v.__get__(None, self)
        return v
    def save(self):
        #print('writing', self.variant_id, 'to database')
        #return Variant.db.variants.update({'variant_id':self.variant_id},self.__dict__,upsert=True)
        pass
    @property
    def kaviar(self):
        if 'kaviar' in self.__dict__: return self.__dict__['kaviar']
    @property
    def status(self):
        return 'M'
    @property
    def HPO(self):
        return []
    @property
    def FILTER(self):
        return self.filter
    @property
    def filter(self):
        self.__dict__['filter']=self.__dict__['FILTER']
        return self.__dict__['filter']
    @property
    def hom_samples(self):
        if 'hom_samples' in self.__dict__: return self.__dict__['hom_samples']
        q=vcf.vcf_query(variant_str=self.variant_id)
        self.__dict__.update(q)
        print(self.save())
        return self.__dict__['hom_samples']
    @property
    def het_samples(self):
        if 'het_samples' in self.__dict__: return self.__dict__['het_samples']
        q=vcf.vcf_query(variant_str=self.variant_id)
        self.__dict__.update(q)
        print(self.save())
        return self.__dict__['het_samples']
    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
    def get_minimal_representation(self): 
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
        pos = int(self.pos)
        # If it's a simple SNV, don't remap anything
        if len(self.ref) == 1 and len(self.alt) == 1: return self.pos, self.ref, self.alt
        # strip off identical suffixes
        while(self.alt[-1] == self.ref[-1] and min(len(self.alt),len(self.ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(self.alt[0] == self.ref[0] and min(len(self.alt),len(self.ref)) > 1):
            alt = self.alt[1:]
            self.ref = self.ref[1:]
            self.pos += 1
        return self.pos, self.ref, self.alt 
    def add_consequence_to_variant(self):
        worst_csq = worst_csq_with_vep(variant['vep_annotations'])
        if worst_csq is None: return
        variant['major_consequence'] = worst_csq['major_consequence']
        variant['HGVSp'] = get_protein_hgvs(worst_csq)
        variant['HGVSc'] = get_transcript_hgvs(worst_csq)
        variant['HGVS'] = get_proper_hgvs(worst_csq)
        variant['CANONICAL'] = worst_csq['CANONICAL']
        variant['flags'] = get_flags_from_variant(variant)
        if csq_order_dict[variant['major_consequence']] <= csq_order_dict["frameshift_variant"]:
            variant['category'] = 'lof_variant'
        elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["missense_variant"]:
            # Should be noted that this grabs inframe deletion, etc.
            variant['category'] = 'missense_variant'
        elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["synonymous_variant"]:
            variant['category'] = 'synonymous_variant'
        else:
            variant['category'] = 'other_variant'
    @property
    def data(self): return self.__dict__
    @property
    def consequence(self):
        """
        Return the most severe consequence
        """
        if 'consequence' in self.__dict__: return self.__dict__['consequence']
        if 'major_consequence' in self.__dict__: return self.__dict__['major_consequence']
        if 'most_severe_consequence' in self.__dict__: return self.__dict__['most_severe_consequence']
        url='http://grch37.rest.ensembl.org/vep/human/hgvs/%s?content-type=application/json' % self.hgvs.replace('chr','')
        r=requests.get(url)
        print(url)
        d=r.json()
        #if not isinstance(d,list) and len(d) < 1: return None
        if 'error' in d: return None
        d=d[0]
        print(d['most_severe_consequence'])
        self.__dict__['consequence']=d['most_severe_consequence']
        print(self.save())
        return self.__dict__['consequence']
    @property
    def transcripts(self):
        if 'transcript_consequences' in self.__dict__: return [x for x in self.transcript_consequences]
        if 'transcripts' in self.__dict__: return self.__dict__['transcripts']
        url='http://grch37.rest.ensembl.org/vep/human/hgvs/%s?content-type=application/json' % self.hgvs.replace('chr','')
        r=requests.get(url)
        print(url)
        d=r.json()
        if 'error' in d: return None
        d=d[0]
        self.__dict__['transcripts']=list(set([csq['transcript_id'] for csq in d.transcript_consequences]))
        self.__dict__['genes']=list(set([csq['gene_id'] for csq in d.transcript_consequences]))
        print(self.save())
        if not isinstance(d,list) and len(d) < 1: return None
        return self.__dict__['transcripts']
    @property
    def genes(self):
        if 'genes' in self.__dict__: return list(set(self.__dict__['genes']))
        url='http://grch37.rest.ensembl.org/vep/human/hgvs/%s?content-type=application/json' % self.hgvs.replace('chr','')
        r=requests.get(url)
        print(url)
        d=r.json()[0]
        self.__dict__['genes']=list(set([csq['gene_id'] for csq in d['transcript_consequences']]))
        print(self.save())
        return self.__dict__['genes']
    @property
    def canonical_hgvsp(self):
        if 'canonical_hgvsp' in self.__dict__:
            return self.__dict__['canonical_hgvsp']
        else:
            return []
    @property
    def canonical_hgvsc(self):
        if 'canonical_hgvsc' in self.__dict__:
            return self.__dict__['canonical_hgvsc']
        else:
            return []
    @property
    def p_hgvs(self):
        """
        Takes consequence dictionary, returns proper variant formatting for synonymous variants
        """
        if '%3D' in csq['HGVSp']:
            try:
                amino_acids = ''.join([protein_letters_1to3[x] for x in csq['Amino_acids']])
                return "p." + amino_acids + csq['Protein_position'] + amino_acids
            except Exception, e:
                print 'Could not create HGVS for: %s' % csq
        return csq['HGVSp'].split(':')[-1]
    @property
    def snpeff(self):
        if 'snpeff' in self.__dict__: return self.__dict__['snpeff']
        self.__dict__['snpeff'] = rest.mv.getvariant('chr%s:g.%s%s>%s'%(self.chrom,self.pos,self.ref,self.alt,),fields='snpeff')
        return self.__dict__['snpeff']
    @property
    def canonical_cadd(self):
        if 'canonical_cadd' in self.__dict__: return self.__dict__['canonical_cadd']
        return ''
    @property
    def cadd(self):
        if 'canonical_cadd' in self.__dict__: return self.__dict__['canonical_cadd']
        if 'cadd' in self.__dict__: return self.__dict__['cadd'].get('phred',None)
        cadd = rest.mv.getvariant('chr%s:g.%s%s>%s'%(self.chrom,self.pos,self.ref,self.alt,),fields='cadd')
        if cadd and 'cadd' in cadd:
            self.__dict__['cadd']=cadd['cadd']
        else:
            self.__dict__['cadd']={}
        print(self.save())
        return self.__dict__['cadd'].get('phred',None)
    @property
    def vep_annotations(self):
        if  'vep_annotations' in self.__dict__: return self.__dict__['vep_annotations']
        if  'transcript_consequences' in self.__dict__: return self.__dict__['transcript_consequences']
        self.__dict__['vep_annotations']=rest.vep_anno(self.chrom, self.pos, self.ref, self.alt)
        print('number of transcripts:', len(self.__dict__['vep_annotations']))
        self.__dict__['transcript_consequences']=self.__dict__['vep_annotations'][0]['transcript_consequences']
        self.__dict__['gene_name_upper']=self.__dict__['transcript_consequences'][0]['gene_symbol']
        print('gene_symbol', self.__dict__['gene_name_upper'])
        #print(self.__dict__['vep_annotations'])
        #self.__dict__['vep_annotations'] = order_vep_by_csq(self.__dict__['vep_annotations']) 
        #self.ordered_csqs = [x['major_consequence'] for x in self.__dict__['vep_annotations']]
        # Close but not quite there
        #ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',')
        #consequences = defaultdict(lambda: defaultdict(list))
        #for annotation in self.data['vep_annotations']:
            #annotation['HGVS'] = get_proper_hgvs(annotation)
            #consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
        return self.__dict__['vep_annotations']
    @property
    def transcript_consequences(self):
        if 'transcript_consequences' in self.__dict__: return self.__dict__['transcript_consequences']
        #print(self.vep_annotations)
        return self.__dict__['transcript_consequences']
    @vep_annotations.setter
    def vep_annotations(self,value):
        self.__dict__['vep_annotations']=value
    @property
    def in_exac(self):
        if 'EXAC' in self.__dict__ and self.__dict__['EXAC'] and len(self.__dict__['EXAC'])>0:
            self.__dict__['in_exac']=True
        else:
            self.__dict__['in_exac']=False
        return self.__dict__['in_exac']
    @property
    def EXAC(self):
        if 'EXAC' in self.__dict__ and self.__dict__['EXAC']:
            if ',' in str(self.__dict__['EXAC']['AC_Hom']):
                self.multi=True
                self.__dict__['EXAC']['AC_Hom']=int(self.__dict__['EXAC']['AC_Hom'].strsplit(',')[0])
                print self.variant_id
            self.__dict__['EXAC']['total_homs']=float(self.__dict__['EXAC']['AC_Hom'])/2
            return self.__dict__['EXAC']
        else:
            return None
        if 'EXAC_freq' in self.__dict__: return self.__dict__['EXAC_freq']
        #self.__dict__['EXAC_freq']=rest.exac_anno(self.data['variant_id'],update=False)
        if len(self.__dict__['EXAC_freq'])>0:
           self.__dict__['in_exac']=True
        else:
           self.__dict__['in_exac']=False
        #print(self.save())
        return self.__dict__['EXAC_freq']
    @EXAC.setter
    def EXAC(self,value):
        self.__dict__['ExAC_freq']=value
        return self.__dict__['EXAC_freq']
    @property
    def ExAC_freq(self):
        if 'ExAC_freq' in self.__dict__ and 'total_homs' in self.__dict__['ExAC_freq']: return self.__dict__['ExAC_freq']
        self.__dict__['ExAC_freq']=rest.exac_anno(self.variant_id,update=False)
        print(self.__dict__['ExAC_freq'].keys())
        #print(self.save())
        return self.__dict__['ExAC_freq']
    @property
    def WT_COUNT(self):
        if 'WT_COUNT' in self.__dict__: return self.__dict__['WT_COUNT']
        q=vcf.vcf_query(variant_str=self.variant_id)
        if q is None: raise Exception('ERROR',self.variant_id)
        self.__dict__.update(q)
        print(self.save())
        return self.__dict__['WT_COUNT']
    @property
    def HOM_COUNT(self):
        if 'HOM_COUNT' in self.__dict__: return self.__dict__['HOM_COUNT']
        q=vcf.vcf_query(variant_str=self.variant_id)
        if q is None: raise Exception('ERROR',self.variant_id)
        self.__dict__.update(q)
        print(self.save())
        return self.__dict__['HOM_COUNT']
    @property
    def allele_num(self):
        if 'allele_num' in self.__dict__: return self.__dict__['allele_num']
        q=vcf.vcf_query(variant_str=self.variant_id)
        if q is None: raise Exception('ERROR',self.variant_id)
        self.__dict__.update(q)
        print(self.save())
        return self.__dict__['allele_num']
    def get_flags_from_variant(self):
        flags = []
        if 'mnps' in variant:
            flags.append('MNP')
        #lof_annotations = [x for x in variant['vep_annotations'] if x['LoF'] != '']
        lof_annotations = []
        if not len(lof_annotations): return flags
        if all([x['LoF'] == 'LC' for x in lof_annotations]):
            flags.append('LC LoF')
        if all([x['LoF_flags'] != '' for x in lof_annotations]):
            flags.append('LoF flag')
        return flags
    @property
    def HUGO(self):
        if 'gene_name_upper' in self.__dict__: return self.__dict__['gene_name_upper']
        if 'canonical_gene_name_upper' in self.__dict__: return self.__dict__['canonical_gene_name_upper'][0]
        else: print(self.variant_id)
        return ''
        #self.vep_annotations
        #print(self.save())
        return self.__dict__['gene_name_upper']
    @property
    def description(self):
        if 'description' in self.__dict__: return self.__dict__['description']
        g=Variant.db.genes.find_one({'gene_name_upper':self.HUGO})
        self.__dict__['description']=g.get('full_gene_name','')
        return self.__dict__['description']
    @property
    def OMIM(self):
        if 'OMIM' in self.__dict__: return self.__dict__['OMIM']
        #self.__dict__['OMIM']=self.vep_annotations[0]['SYMBOL']
        #print(self.save())
        #return self.__dict__['OMIM']
        return ''
    @property
    def p_change(self):
        if 'p_change' in self.__dict__: return self.__dict__['p_change']
        if 'HGVSp' in self.__dict__: return self.__dict__['HGVSp']
        #if 'canonical_hgvsp' in self__dict__: return self.__dict__['canonical_hgvsp']
        self.__dict__['p_change']=dict()
        #self.__dict__['p_change']=
        #trans['hgvsp'].split(':')[1]
        self.__dict__['p_change']['exon']=''
        self.__dict__['p_change']['gene_id']=self.genes[0]
        self.__dict__['p_change']['transcript_id']=self.canonical_transcript[0]
        self.__dict__['p_change']['hgvs_c']=self.canonical_hgvsc[0]
        self.__dict__['p_change']['hgvs_p']=self.canonical_hgvsp[0]
        return self.__dict__['p_change']
    # get db
    def stuff():
        if 'consequence' in self.__dict__ and len(self.__dict__['consequence']): return self.__dict__['consequence']
        pp = pprint.PrettyPrinter(indent=10)
        v['Consequence']=[transcript['consequence_terms'][0] for transcript in v['vep_annotations']['transcript_consequences']]
        v['vep_annotations']['Consequence']=[csq for csq in v['Consequence']]
        print ('CSQ')
        print( v['vep_annotations']['Consequence'] )
        worst_csq = worst_csq_with_vep(variant['vep_annotations'])
        if worst_csq is None: return
        variant['major_consequence'] = worst_csq['major_consequence']
        variant['HGVSp'] = get_protein_hgvs(worst_csq)
        variant['HGVSc'] = get_transcript_hgvs(worst_csq)
        variant['HGVS'] = get_proper_hgvs(worst_csq)
        variant['CANONICAL'] = worst_csq['CANONICAL']
        variant['flags'] = get_flags_from_variant(variant)
        if csq_order_dict[variant['major_consequence']] <= csq_order_dict["frameshift_variant"]:
            variant['category'] = 'lof_variant'
        elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["missense_variant"]:
            # Should be noted that this grabs inframe deletion, etc.
            variant['category'] = 'missense_variant'
        elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["synonymous_variant"]:
            variant['category'] = 'synonymous_variant'
        else:
            variant['category'] = 'other_variant'
    def worst_csq_with_vep(self, annotation_list):
        """
        Takes list of VEP annotations [{'Consequence': 'frameshift', Feature: 'ENST'}, ...]
        Returns most severe annotation (as full VEP annotation [{'Consequence': 'frameshift', Feature: 'ENST'}])
        Also tacks on worst consequence for that annotation (i.e. worst_csq_from_csq)
        :param annotation_list:
        :return worst_annotation:
        """
        if len(annotation_list) == 0: return None
        worst = annotation_list[0]
        for annotation in annotation_list:
            if compare_two_consequences(annotation['Consequence'], worst['Consequence']) < 0:
                worst = annotation
            elif compare_two_consequences(annotation['Consequence'], worst['Consequence']) == 0 and annotation['CANONICAL'] == 'YES':
                worst = annotation
        worst['major_consequence'] = worst_csq_from_csq(worst['Consequence'])
        return worst
    



