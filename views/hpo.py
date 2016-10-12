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


@app.route('/hpo')
def hpo_main():
    # HPO summary page
    # major groups, borrowed from phenotips
    major_groups = {'GROWTH PARAMETERS':['HP:0000256','HP:0000252','HP:0000098','HP:0004322','HP:0004324','HP:0004325','HP:0001508','HP:0001528'],'CRANIOFACIAL':['HP:0001363','HP:0000204','HP:0000175','HP:0001999'],'EYE DEFECTS':['HP:0000505','HP:0000481','HP:0000589','HP:0000593','HP:0000518','HP:0000479','HP:0000587','HP:0000568','HP:0000639','HP:0000486','HP:0000601','HP:0000316'],'EAR DEFECTS':['HP:0000407','HP:0000405','HP:0004467','HP:0000384','HP:0000356','HP:0000359'],'CUTANEOUS':['HP:0000953','HP:0001010','HP:0005306','HP:0011276'],'CARDIOVASCULAR':['HP:0001631','HP:0001629','HP:0001674','HP:0001680','HP:0001636','HP:0001638','HP:0011675'],'RESPIRATORY':['HP:0000776','HP:0002088'],'MUSCULOSKELETAL':['HP:0002652','HP:0002659','HP:0009816','HP:0009824','HP:0100490','HP:0001836','HP:0006101','HP:0001770','HP:0100258','HP:0100259','HP:0001180','HP:0001849','HP:0002650','HP:0000925','HP:0001371','HP:0001762'],'GASTROINTESTINAL':['HP:0002032','HP:0002575','HP:0001543','HP:0001539','HP:0002251','HP:0001396','HP:0002910','HP:0001738','HP:0000819'],'GENITOURINARY':['HP:0000107','HP:0000085','HP:0000069','HP:0000795','HP:0000062','HP:0000047','HP:0000028'],'BEHAVIOR, COGNITION AND DEVELOPMENT':['HP:0001263','HP:0010862','HP:0002194','HP:0000750','HP:0001328','HP:0001256','HP:0002342','HP:0010864','HP:0007018','HP:0000717','HP:0000708'],'NEUROLOGICAL':['HP:0001290','HP:0001250','HP:0001251','HP:0001332','HP:0002072','HP:0001257','HP:0010301','HP:0002011']}
    hpo_freq = lookups.get_hpo_size_freq('hpo_freq.tsv')
    return str(hpo_freq)


@app.route('/hpo/<hpo_id>')
@requires_auth
def hpo_page(hpo_id):
    patients_db=get_db('patients')
    db=get_db()
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    print(hpo_id)
    if not hpo_id.startswith('HP:'):
        hpo_id=hpo_db.hpo.find_one({'name':hpo_id})['id'][0]
    print(hpo_id)
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
    lit_genes=[lookups.get_gene_by_name(db, r['Gene-Name']) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)
    print('num patients', len(patients))
    pmids=[]
    obs_genes=[]
    if len(patients) < 100:
        for p in patients:
                p2=db.patients.find_one({'external_id':p['external_id']},{'rare_variants':1})
                if not p2: continue
                #print(p['external_id'])
                for rv in p2.get('rare_variants',[]):
                    if 'HUGO' in rv: obs_genes+=[rv['HUGO']]
    obs_genes=Counter(obs_genes)
    obs_genes=obs_genes.most_common(10)
    #obs_genes=[g for g in obs_genes if obs_genes[g]>5]
    obs_genes=[g for g in obs_genes]
    # candidate genes
    candidate_genes = [p.get('genes',[]) for p in patients]
    # solved genes
    solved_genes = [p.get('solved',[]) for p in patients]
    #print(intersect(obs_genes.keys(),lit_genes))
    #print(Counter([rv['HUGO'] for rv in db.patients.find_one({'external_id':p['external_id']},{'rare_variants':1})]['rare_variants']))
    ## only return common variants if there are many individuals
    ##rsession.voidEval('common_variants <- common.variants')
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
    #vcf_reader = pysam.VariantFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % '22')
    #for record in vcf_reader:
        #for s in external_ids:
            #r=record.samples[s]
            #if 'GT' in r: print(r['GT'])
    return render_template('hpo.html',hpo_id=hpo_id,hpo_name=hpo_name,
            individuals=[str(p['external_id']) for p in patients],
            lit_genes=lit_genes,
            obs_genes=obs_genes,
            hpo_gene=hpo_gene,
            pmids=pmids,variants=[])


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


