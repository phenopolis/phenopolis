#! /bin/env python
from __future__ import print_function
import sys
import pymongo
from optparse import OptionParser


usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("--host", dest="host", help="", default='localhost')
parser.add_option("--port", dest="port", help="port", default=27017)
parser.add_option("--inheritance-mode", dest="inheritance_mode", help="inheritance mode: dominant or recessive",default='recessive')
parser.add_option("--hpo", dest="hpo", help="")
parser.add_option("--cache", dest="cache", help="cache", default=False, action="store_true")
(options, args) = parser.parse_args()

conn = pymongo.MongoClient(host=options.host, port=options.port)
hpo_id=options.hpo


def phenogenon(hpo_id,lit_genes,omim_genes,recessive_genes,dominant_genes,cache=False):
    cache_db=conn['cache']
    temp=cache_db.phenogenon_cache.find_one({'hpo_id':hpo_id})
    if temp and cache:
        lit_genes.extend(temp['lit_genes'])
        omim_genes.extend(temp['omim_genes'])
        recessive_genes.extend(temp['recessive_genes'])
        dominant_genes.extend(temp['dominant_genes'])
        return
    hpo_db=conn['hpo']
    db=conn['uclex']
    def f(r):
        g=db.genes.find_one({'gene_name_upper':r['Gene-Name'].upper()},{'_id':0})
        if not g: return
        phenogenon=db.gene_hpo.find_one({'gene_id':g['gene_id']})
        if not phenogenon: return
        het=phenogenon.get('het',{}).get(hpo_id,{})
        hom_comp=phenogenon.get('hom_comp',{}).get(hpo_id,{})
        if 'data' in het: del het['data']
        if 'data' in hom_comp: del hom_comp['data']
        g['phenogenon']={ 'het':het, 'hom_comp': hom_comp}
        return g
    lit_genes=[f(r) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    lit_genes=[lg for lg in lit_genes if lg]
    omim_genes.extend(map(lambda x: x['gene_id'], lit_genes))
    phenogenon=db.hpo_gene.find_one({'hpo_id':hpo_id})
    if phenogenon: phenogenon=phenogenon['data']['unrelated']
    else: phenogenon={'recessive':[],'dominant':[]}
    recessive_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'],'known':x['gene_id'] in omim_genes} for x in phenogenon['recessive']])
    dominant_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'], 'known':x['gene_id'] in omim_genes} for x in phenogenon['dominant']])
    #print({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})
    cache_db.phenogenon_cache.insert_one({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})


def phenogenon2(hpo_id,lit_genes,omim_genes,recessive_genes,dominant_genes,cache=False):
    cache_db=conn['cache']
    temp=cache_db.phenogenon_cache.find_one({'hpo_id':hpo_id})
    if temp and cache:
        lit_genes.extend(temp['lit_genes'])
        omim_genes.extend(temp['omim_genes'])
        recessive_genes.extend(temp['recessive_genes'])
        dominant_genes.extend(temp['dominant_genes'])
        return
    hpo_db=conn['hpo']
    db=conn['uclex']
    def f(r):
        g=db.genes.find_one({'gene_name_upper':r['Gene-Name'].upper()},{'_id':0})
        if not g: return
        phenogenon=db.gene_hpo.find_one({'gene_id':g['gene_id']})
        if not phenogenon: return
        het=phenogenon.get('het',{}).get(hpo_id,{})
        hom_comp=phenogenon.get('hom_comp',{}).get(hpo_id,{})
        if 'data' in het: del het['data']
        if 'data' in hom_comp: del hom_comp['data']
        g['phenogenon']={ 'het':het, 'hom_comp': hom_comp}
        return g
    lit_genes=[f(r) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    lit_genes=[lg for lg in lit_genes if lg]
    omim_genes.extend(map(lambda x: x['gene_id'], lit_genes))
    phenogenon=db.hpo_gene.find_one({'hpo_id':hpo_id})
    if phenogenon: phenogenon=phenogenon['data']['unrelated']
    else: phenogenon={'recessive':[],'dominant':[]}
    recessive_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'],'known':x['gene_id'] in omim_genes} for x in phenogenon['recessive']])
    dominant_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'], 'known':x['gene_id'] in omim_genes} for x in phenogenon['dominant']])
    #print({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})
    cache_db.phenogenon_cache.insert_one({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})



lit_genes=[]
omim_genes=[]
recessive_genes=[]
dominant_genes=[]
if options.inheritance_mode=='dominant':
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes,options.cache)
    if len(dominant_genes)==0: sys.exit(0)
    print(','.join([k for k in dominant_genes[0].keys()]))
    for g in dominant_genes:
        print(','.join([str(g[k]) for k in g]))
elif options.inheritance_mode=='recessive':
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes,options.cache)
    if len(recessive_genes)==0: sys.exit(0)
    print(','.join([k for k in recessive_genes[0].keys()]))
    for g in recessive_genes:
        print(','.join([str(g[k]) for k in g]))
elif options.inheritance_mode=='literature':
    phenogenon2(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes, options.cache)
    names=['gene_name','phenogenon.dominant_pvalue','phenogenon.recessive_pvalue']
    print(','.join([k for k in names]))
    for g in lit_genes:
        gene_name=g['gene_name']
        dominant_pvalue=str(g['phenogenon']['het'].get('unrelated_dominant_all_p_val',None))
        recessive_pvalue=str(g['phenogenon']['hom_comp'].get('unrelated_recessive_p_val',None))
        print(','.join([gene_name,dominant_pvalue,recessive_pvalue]))
        
        
