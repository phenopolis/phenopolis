from views import *
from lookups import *
import rest as annotation
import requests
import primer3
import myvariant
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import phizz
import random
import orm
import vcf


def phenogenon(hpo_id,lit_genes,omim_genes,recessive_genes,dominant_genes):
    hpo_db=get_db('hpo')
    db=get_db()
    for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id}):
        g=db.genes.find_one({'gene_name_upper':r['Gene-Name'].upper()},{'_id':0})
        if not g: continue
        phenogenon=db.gene_hpo.find_one({'gene_id':g['gene_id']})
        if not phenogenon: continue
        g['phenogenon']={
                'het':phenogenon.get('het',{}).get(hpo_id,{}),
                'hom_comp':phenogenon.get('hom_comp',{}).get(hpo_id,{})
                }
        lit_genes+=[g]
    omim_genes.extend(map(lambda x: x['gene_id'], lit_genes))
    phenogenon=db.hpo_gene.find_one({'hpo_id':hpo_id})
    if phenogenon: phenogenon=phenogenon['data']['unrelated']
    else: phenogenon={'recessive':[],'dominant':[]}
    recessive_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'],'known':x['gene_id'] in omim_genes} for x in phenogenon['recessive']])
    dominant_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'], 'known':x['gene_id'] in omim_genes} for x in phenogenon['dominant']])


@app.route('/phenogenon_json/<hpo_id>')
@requires_auth
def phenogenon_json(hpo_id):
    print 'PHENOGENON_JSON'
    #print(intersect(obs_genes.keys(),lit_genes))
    #print(Counter([rv['HUGO'] for rv in db.patients.find_one({'external_id':p['external_id']},{'rare_variants':1})]['rare_variants']))
    ## only return common variants if there are many individuals
    ##rsession.voidEval('common_variants <- common.variants')
    lit_genes=[]
    omim_genes=[]
    recessive_genes=[]
    dominant_genes=[]
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes)
    print(len(lit_genes))
    print(len(omim_genes))
    print(len(recessive_genes))
    print(len(dominant_genes))
    return jsonify( result={
                    'lit_genes':lit_genes,
                    'omim_genes':omim_genes,
                    'recessive_genes':recessive_genes,
                    'dominant_genes':dominant_genes} )

@app.route('/hpo/<hpo_id>')
@requires_auth
def hpo_page(hpo_id):
    patients_db=get_db('patients')
    db=get_db()
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    print hpo_id 
    if not hpo_id.startswith('HP:'):
        hpo_term=hpo_db.hpo.find_one({'name':re.compile('^'+hpo_id+'$',re.IGNORECASE)})
        if not hpo_term: return hpo_id+' does not exist'
        hpo_id=hpo_term['id'][0]
    print hpo_id 
    hpo_name=hpo_db.hpo.find_one({'id':hpo_id})['name'][0]
    print('HPO ANCESTORS')
    hpo_ancestors=lookups.get_hpo_ancestors(hpo_db,hpo_id)
    hpo_gene = db.hpo_gene.find_one({'hpo_id':hpo_id})
    print(len(hpo_ancestors))
    print([h['name'] for h in hpo_ancestors])
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    #r=patients_db.hpo.find_one({'hp_id':hpo_id})
    #if r: external_ids=r['external_ids']
    #else: external_ids=[]
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)
    print('num patients', len(patients))
    pmids=[]
    obs_genes=[]
    #if len(patients) < 100:
    if False:
        for p in patients:
                p2=db.patients.find_one({'external_id':p['external_id']},{'rare_variants':1})
                if not p2: continue
                for rv in p2.get('rare_variants',[]):
                    if 'HUGO' in rv: obs_genes+=[rv['HUGO']]
    obs_genes=Counter(obs_genes)
    obs_genes=obs_genes.most_common(10)
    print(obs_genes)
    #obs_genes=[g for g in obs_genes if obs_genes[g]>5]
    obs_genes=[g for g in obs_genes]
    # candidate genes
    candidate_genes = [p.get('genes',[]) for p in patients]
    # solved genes
    solved_genes = [p.get('solved',[]) for p in patients]
    ## private variants (not seen in others in the cohort)
    ##rsession.voidEval('common_variants <- common.variants')
    #variants=rsession.r.private_variants(hpo_patients)
    #if type(variants) is str:
        #variants=[variants]
    #else:
        #variants=variants.tolist()
    #print('num variants',len(variants),)
    #variants=[db.variants.find_one({'variant_id':v.replace('_','-')}) for v in variants[:100]]
    #[variant for variant in lookups.get_variants_in_gene(db, g['gene_id'])]
       #if variant['major_consequence']!='stop_gained': continue
       #print(variant)
       #break
    #print( lookups.get_variants_in_gene(db, 'CNNM4') )
    #vcf_reader = pysam.VariantFile('chr%s.vcf.gz' % '22')
    #for record in vcf_reader:
        #for s in external_ids:
            #r=record.samples[s]
            #if 'GT' in r: print(r['GT'])
    hpo_db=get_db('hpo')
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=db.patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        solved_patient=db.solved_patients.find_one({'external_id':p['external_id']})
        if solved_patient and session['user']!='demo': p['solved_variants']=solved_patient.get('genes',{})
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    #eids=[p['eid'] for p in patients]
    #print(eids)
    #patients=get_db('patients').patients.find({'external_id':{'$in':eids}})
    #patients=get_db('patients').patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    patients=[f(p) for p in patients[:500] if 'external_id' in p]
    #print recessive_genes
    #print dominant_genes
    lit_genes=[r['Gene-Name'] for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    skat_genes=db.skat.find({'HPO':hpo_id})
    skat_genes=[g for g in skat_genes if g['FisherPvalue']<0.05 and g['SKATO']<0.005]
    for g in skat_genes:
        pli=get_db('exac').pli.find_one({'gene':g['Symbol']})
        if pli:
            g['pli']=pli['pLI'] 
        else:
            g['pli']=-1
    return render_template('hpo.html',
            title=hpo_id,
            hpo_id=hpo_id,
            hpo_name=hpo_name,
            individuals=[p for p in patients],
            lit_genes=lit_genes,
            obs_genes=obs_genes,
            recessive_genes=[],
            dominant_genes=[],
            hpo_gene=hpo_gene,
            skat_genes=skat_genes,
            pmids=pmids,
            variants=[])


@auth.verify_password
def verify_pw(username,password): return check_auth(username, password)

@app.route('/hpo_json/<hpo_id>')
@auth.login_required
def hpo_json(hpo_id):
    patients_db=get_db('patients')
    db=get_db()
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    print(hpo_id)
    if not hpo_id.startswith('HP:'):
        hpo_term=hpo_db.hpo.find_one({'name':hpo_id})
        hpo_id=hpo_term['id'][0]
    print(hpo_id)
    hpo_term=hpo_db.hpo.find_one({'id':hpo_id})
    hpo_name=hpo_term['name'][0]
    if 'is_a' in hpo_term:
        parents=[ pid for pid in hpo_term['is_a'] ]
    else:
        parents=[]
    print('HPO ANCESTORS')
    hpo_ancestors=lookups.get_hpo_ancestors(hpo_db,hpo_id)
    #print(lookups.get_hpo_ancestors_array(hpo_db,hpo_id))
    print(hpo_ancestors)
    print(len(hpo_ancestors))
    print([h['name'] for h in hpo_ancestors])
    #hpo_ancestors=dict((h['id'][0],h['name'][0]) for h in hpo_ancestors)
    hpo_ancestors=[{'hpo_id':h['id'][0],'hpo_name':h['name'][0]} for h in hpo_ancestors]
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    #r=patients_db.hpo.find_one({'hp_id':hpo_id})
    #if r: external_ids=r['external_ids']
    #else: external_ids=[]
    genes=[lookups.get_gene_by_name(db, r['Gene-Name']) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    print('num genes', len(genes))
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)
    print('num patients', len(patients))
    #return jsonify(result={'hpo_id':hpo_id,'hpo_name':hpo_name,'individuals':[str(p['external_id']) for p in patients],'genes':genes})
    return jsonify(result={'hpo_id':hpo_id,'hpo_name':hpo_name,'individuals':[str(p['external_id']) for p in patients],'hpo_ancestors':hpo_ancestors, 'parents':parents})


