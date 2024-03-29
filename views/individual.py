from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
from config import config
if config.IMPORT_PYSAM_PRIMER3:
    import pysam
import csv
#hpo lookup
import orm
from pprint import pprint
import os
import json
import pymongo
import sys
import re
import itertools
from urllib2 import HTTPError, URLError
import csv
from collections import defaultdict, Counter
#import rest as annotation
from optparse import OptionParser
import mygene
import lookups
from orm import Patient
import requests

from neo4j.v1 import GraphDatabase, basic_auth

def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
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
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


@app.route('/update_patient_data/<individual>',methods=['POST'])
@requires_auth
def update_patient_data(individual):
    if session['user']=='demo': return 'not permitted'
    print(request.form)
    consanguinity=request.form.getlist('consanguinity_edit[]')[0]
    gender=request.form.getlist('gender_edit[]')[0]
    genes=request.form.getlist('genes[]')
    features=request.form.getlist('feature[]')
    print('INDIVIDUAL',individual)
    print('GENDER',gender)
    print('CONSANGUINITY',consanguinity)
    print('GENES',genes)
    print('FEATURES',features)
    print(individual)
    external_id=individual
    individual=get_db(app.config['DB_NAME_PATIENTS']).patients.find_one({'external_id':external_id})
    print('edit patient gender')
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'sex':{'female':'F','male':'M','unknown':'U'}[gender]}}))
    print('edit patient genes')
    individual['genes']=[]
    for g in genes:
        gene=get_db(app.config['DB_NAME']).genes.find_one({'gene_name_upper':g})
        print(gene)
        if gene in [g['gene'] for g in individual['genes']]: continue
        if not gene: continue
        individual['genes'].append({'gene':g, 'status':'candidate'})
    print(individual['genes'])
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'genes':individual['genes']}}))
    print('edit patient features')
    individual['features']=[]
    for f in features:
        hpo=get_db(app.config['DB_NAME_HPO']).hpo.find_one({'name':re.compile('^'+f+'$',re.IGNORECASE)})
        if not hpo: continue
        if hpo in [h['label'] for h in individual['features']]: continue
        individual['features'].append({'id':hpo['id'][0], 'label':hpo['name'][0], 'observed':'yes'})
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'features':individual['features']}}))
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'observed_features':[f for f in individual['features'] if f['observed']=='yes']}}))
    print('edit patient consanguinity')
    individual['family_history']=individual.get('family_history',{})
    if (consanguinity)=='unknown':
        individual['family_history']['consanguinity']=None
    elif consanguinity.lower()=='yes':
        individual['family_history']['consanguinity']=True
    elif consanguinity.lower()=='no':
        individual['family_history']['consanguinity']=False
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'family_history':individual['family_history']}}))
    # also trigger refresh of that individual for individuals summary page
    patient=Patient(external_id,get_db(app.config['DB_NAME_PATIENTS']))
    print(patient.consanguinity)
    print(patient.observed_features)
    print(patient.genes)
    print(patient.gender)
    individuals_update([external_id])
    return jsonify({'success': True}), 200


@app.route('/individual_json/<individual>')
@requires_auth
def individual_json(individual):
    patient=Patient(individual,patient_db=get_db(app.config['DB_NAME_PATIENTS']))
    #PP.addPatientGenderInfo(data.result.sex); 
    #PP.addPatientFeaturesInfo(data.result.observed_features);
    #PP.addPatientConsanguinityInfo(data.result.family_history);
    #PP.addPatientGenesInfo(data.result.genes);
    #PP.submitEditedIndividual(patientId);
    return patient.json()


@app.route('/individual/<individual>')
@requires_auth
#@cache.cached(timeout=24*3600)
def individual_page(individual):
    patient=Patient(individual,patient_db=get_db(app.config['DB_NAME_PATIENTS']),variant_db=get_db(app.config['DB_NAME']),hpo_db=get_db(app.config['DB_NAME_HPO']))
    #if session['user']=='demo': individual=decrypt(str(individual))
    # make sure that individual is accessible by user
    if not lookup_patient(db=get_db(app.config['DB_NAME_USERS']),user=session['user'],external_id=individual): return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    # TODO
    # mode of inheritance in hpo terms: HP:0000005
    #print lookups.get_hpo_children(hpo_db, 'HP:0000005')
    #patient['global_mode_of_inheritance']=patient2.get('global_mode_of_inheritance',None)
    # minimise it
    patient.__dict__['hpo_ids']=lookups.hpo_minimum_set(get_db(app.config['DB_NAME_HPO']), patient.hpo_ids)
    hpo_gene=get_hpo_gene(patient.hpo_ids)
    # get pubmedbatch scores
    pubmedbatch = {}
    genes = {}
    # is this still updating?
    if type(pubmedbatch) is dict:
        update_status = pubmedbatch.get('status', 0)
    else:
        update_status=0
    # get known and retnet genes
    known_genes=[x['gene_name'] for x in db.retnet.find()]
    RETNET = dict([(i['gene_name'],i) for i in db.retnet.find({},projection={'_id':False})])
    print 'get pubmed score and RETNET'
    gene_info=dict()
    individuals=dict()
    #
    genes=[]
    #genes['homozygous_variants']=[v['canonical_gene_name_upper'] for v in patient.homozygous_variants]
    #genes['compound_hets']=[v['canonical_gene_name_upper'] for v in patient.compound_het_variants]
    #genes['rare_variants']=[v['canonical_gene_name_upper'] for v in patient.rare_variants]
            # print(g, genes_pubmed[g])
    # figure out the order of columns from the variant row
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?.*>",file('templates/individual-page-tabs/individual_variant_row.tmpl','r').read())
    if session['user']=='demo': table_headers=table_headers[:-1]
    print table_headers
    # get a list of genes related to retinal dystrophy. only relevant to subset group of ppl. talk to Jing or Niko for other cohorts. Note that dominant p value only counts paitents with 1 qualified variant on the gene. 
    # current setting: unrelated, exac_af 0.01 for recessive, 0.001 for dominant, cadd_phred 15
    print 'get phenogenon genes'
    retinal_genes = {}
    return render_template('individual.html', 
            patient=patient,
            table_headers=table_headers,
            pubmedbatch=pubmedbatch,
            pubmed_db=get_db('pubmed_cache'),
            genes = genes,
            individuals=individuals,
            hpo_gene = hpo_gene,
            gene_info={},
            update_status = 0,
            retinal_genes = {},
            feature_venn = [])


def get_feature_venn(patient):
    s="""
    MATCH (p:Person)-[:PersonToObservedTerm]->(t:Term)--(g:Gene)
    WHERE p.personId='%s'
    RETURN t.termId, t.name, g.gene_id, g.gene_name
    """ % patient
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)

    data = []
    for r in result:
        data.append({
            'hpo_id': r['t.termId'],
            'hpo_term': r['t.name'],
            'gene_id': r['g.gene_id'],
            'gene_name': r['g.gene_name']
        })

    hpo_terms=[(k,v,) for k, v, in dict([(x['hpo_id'],x['hpo_term'],) for x in data]).items()]
    hpo_gene=dict()
    for x in data:
        hpo_gene[x['hpo_id']]=hpo_gene.get(x['hpo_id'],[])+[x['gene_name']]

    genes = {}
    feature_combo = []
    feature_venn = []
    print "get combinatorics of features to draw venn diagram"
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = [hpo_terms[i][1] for i in combo]
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':list(frozenset(hpo_gene.get(x,"")))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = list(frozenset(feature_venn[-1]['value']) & frozenset(hpo_gene[hpo_terms[combo[ind]][0]]))
    return feature_venn



@app.route('/venn_json/<individual>')
@requires_auth
def venn_json(individual):
    feature_venn=get_feature_venn(individual)
    return jsonify(result=feature_venn)


def patient_variants():
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_het_variants']:
        for v in patient.__dict__[mm]:
            if 'canonical_gene_name_upper' not in v: v['canonical_gene_name_upper']=v['Gene']
            gene=v['canonical_gene_name_upper']
            pubmed_key = '_'.join([gene,patient.get('pubmed_key','')])
            gene_info[gene]=dict()
            if gene in known_genes: 
                gene_info[gene]['known']=True
                pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if gene not in RETNET: continue
            gene_info[gene]['disease'] = RETNET[gene]['disease']
            gene_info[gene]['omim'] = RETNET[gene]['omim']
            gene_info[gene]['mode'] = RETNET[gene]['mode']
            pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if mm != 'rare_variants' or ('d' in gene_info[gene]['mode'] and mm == 'rare_variants') :
                pubmedbatch[pubmed_key] = max(100,pubmedbatch[pubmed_key])
                if gene=='DRAM2':
                    print pubmed_key
                    print pubmedbatch[pubmed_key]
            if 'het_samples' not in v: print(v)
            for s in v['het_samples']:
                if v['HET_COUNT'] < 10:
                    individuals[s]=individuals.get(s,[])+[v]



def get_hpo_gene(hpo_ids):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    hpo_terms = [(i, hpo_db.hpo.find_one({'id':i})['name'][0]) for i in hpo_ids]
    # this has missing HPO ids. see IRDC_batch2_OXF_3001 and #HP:0000593
    hpo_gene=dict()
    for hpo_id,hpo_term, in hpo_terms:
        hpo_gene[hpo_id] = []
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            #gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
            hpo_gene[hpo_id]=hpo_gene.get(hpo_id,[])+[gene_name]
    for k in hpo_gene: hpo_gene[k]=list(frozenset(list(hpo_gene[k])))
    return hpo_gene


def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item

def exomiser(individual):
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, individual, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    x['exomiser']=[]
    for g in list(set(x['genes'])):
        r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
        if not r or not r['EntrezGene ID']: continue
        x['entrezgeneid']=r['EntrezGene ID']
        #url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        url='http://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        if isinstance(r.json(),list):
            x['exomiser']+=r.json()[0]['results']
        else:
            x['exomiser']+=r.json()['results']
    if len(x['exomiser'])<1: x['exomiser']=[{'score':-1}]
    exomiser_scores=[xx['score'] for xx in x['exomiser']]
    i=exomiser_scores.index(max(exomiser_scores))
    x['exomiser']=x['exomiser'][i]


@app.route('/homozygous_variants_json/<individual>')
@requires_auth
def homozgous_variants(individual):
    patient=Patient(individual,patient_db=get_db(app.config['DB_NAME_PATIENTS']),variant_db=get_db(app.config['DB_NAME']),hpo_db=get_db(app.config['DB_NAME_HPO']))
    return jsonify(result=patient.homozygous_variants)


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

@app.route('/homozygous_variants_json2/<individual>')
@requires_auth
def homozygous_variants2(individual):
    allele_freq=float(request.args.get('allele_freq',0.001))
    kaviar_AF=float(request.args.get('kaviar_AF',0.001))
    s=""" MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (t)--(g:Gene)--(gv:GeneticVariant)-[:HomVariantToPerson]-(p:Person), 
    (gv)--(tv:TranscriptVariant)
    WHERE p.personId='%s' AND gv.kaviar_AF < %f AND gv.allele_freq < %f
    WITH gv, g, t, tv
    OPTIONAL
    MATCH
    (gv)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (gv)-[:HomVariantToPerson]-(p3:Person)
    RETURN gv,
    collect(distinct g),
    collect(distinct t),
    collect(distinct tv),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq,)
    print(s)
    with neo4j_driver.session() as neo4j_session: 
        result=neo4j_session.run(s)
    return jsonify(result=[merge_dicts(dict(r[0]),
        {'genes':[dict(x) for x in r[1]]},
        {'terms':[dict(x) for x in r[2]]},
        {'transcript_variants':[dict(x) for x in r[3]]},
        {'het_individuals':[dict(x) for x in r[4]]},
        {'hom_individuals':[dict(x) for x in r[5]]}
        ) for r in result])


@app.route('/compound_het_variants_json2/<individual>',methods=['GET','POST'])
@requires_auth
def compound_het_variants2(individual):
    kaviar_AF=float(request.args.get('kaviar_AF',0.01))
    allele_freq=float(request.args.get('allele_freq',0.01))
    s="""
    MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (g:Gene)--(gv:GeneticVariant)-[:HetVariantToPerson]-(p:Person)
    WHERE p.personId='%s' AND gv.kaviar_AF<%f and gv.allele_freq < %f
    WITH g, collect(distinct gv) AS cgv
    WHERE length(cgv) > 1
    UNWIND cgv as v
    OPTIONAL
    MATCH
    (v)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (v)-[:HomVariantToPerson]-(p3:Person)
    RETURN v,
    collect(distinct g),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq)
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0]),
        {'terms':[]},
        {'genes':[dict(x) for x in r[1]]},
        {'transcript_variants':[]},
        {'het_individuals':[dict(x) for x in r[2]]},
        {'hom_individuals':[dict(x) for x in r[3]]}
        ) for r in result])
 
@app.route('/compound_het_variants_json/<individual>')
@requires_auth
def compound_het_variants(individual):
    patient=Patient(individual,patient_db=get_db(app.config['DB_NAME_PATIENTS']),variant_db=get_db(app.config['DB_NAME']),hpo_db=get_db(app.config['DB_NAME_HPO']))
    return jsonify(result=patient.compound_het_variants)

@app.route('/rare_variants_json2/<individual>')
@requires_auth
def rare_variants2(individual):
    kaviar_AF=float(request.args.get('kaviar_AF',0.01))
    allele_freq=float(request.args.get('allele_freq',0.01))
    s=""" MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (t)--(g:Gene)--(gv:GeneticVariant)-[:HetVariantToPerson]-(p:Person), 
    (gv)--(tv:TranscriptVariant)
    WHERE p.personId='%s' AND gv.kaviar_AF < %f AND gv.allele_freq < %f
    WITH gv, g, t, tv
    OPTIONAL
    MATCH
    (gv)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (gv)-[:HomVariantToPerson]-(p3:Person)
    RETURN gv,
    collect(distinct g),
    collect(distinct t),
    collect(distinct tv),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq,)
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0]),
        {'genes':[dict(x) for x in r[1]]},
        {'terms':[dict(x) for x in r[2]]},
        {'transcript_variants':[dict(x) for x in r[3]]},
        {'het_individuals':[dict(x) for x in r[4]]},
        {'hom_individuals':[dict(x) for x in r[5]]}
        ) for r in result])
 


def load_patient(individual,auth,pubmed_key,hpo='HP:0000001'):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    db = get_db()
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient_id=individual
    patient={u'features': {u'observed': u'yes', u'type': u'phenotype', u'id': hpo}, 'clinicalStatus': {u'clinicalStatus': u'affected'}, u'ethnicity': {u'maternal_ethnicity': [], u'paternal_ethnicity': []}, u'family_history': {}, u'disorders': [], u'life_status': u'alive', u'reporter': u'', u'genes': [], u'prenatal_perinatal_phenotype': {u'prenatal_phenotype': [], u'negative_prenatal_phenotype': []}, u'prenatal_perinatal_history': {u'twinNumber': u''}, u'sex': u'U', u'solved': {u'status': u'unsolved'}}
    eid=patient_id
    if p: patient.update(p)
    #patient_hpo_terms=','.join([f['id'] for f in patient['features'] if f['observed']=='yes'])
    gene_counter=Counter([var['canonical_gene_name_upper'] for var in patient.rare_variants])
    for var in patient['rare_variants']: var['gene_count']=gene_counter[var['canonical_gene_name_upper']]
    patient["pubmedbatch_status"]=0
    pubmed_key="blindness-macula-macular-pigmentosa-retina-retinal-retinitis-stargardt"
    patient["pubmed_key"]=pubmed_key
    #db.patients.update({'external_id':patient_id}, patient, upsert=True)



@app.route('/individual_update/<individual>')
@requires_auth
def individual_update(individual):
    print 'UPDATE'
    print p
    print get_db(app.config['DB_NAME_PATIENTS']).patients.update({'external_id':individual},{'$set':p})
    print 'DB'
    print get_db(app.config['DB_NAME_PATIENTS']).patients.find_one({'external_id':individual})
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
        return redirect(referrer+'/individual/'+individual)
    else:
        return 'done'


'''
progress bar query
'''
@app.route('/pubmedbatch_progress_bar/<id>')
def pubmedbatch_progress(id):
    user = session.get('user') or app.config['DEFAULT_USER']
    progress_id = user + id
    return jsonify(PROGRESS_BAR[progress_id])

'''
get pubmedbatch cache results based on pubmedkey
'''
@app.route('/pubmedbatch-cache/<pubmedkey>')
def pubmedbatch_getcache(pubmedkey):
    db = get_db('pubmedbatch') 
    result = db.cache.find_one({'key':pubmedkey},{'_id':False})
    if result: return jsonify(result)
    else: return jsonify('')


@app.route('/homozygous_individuals_json/<variant_id>')
@requires_auth
def get_homozygous_individuals(variant_id):
    s="""
    MATCH
    (v)-[:HomVariantToPerson]-(p:Person)
    WHERE v.variantId='%s'
    RETURN p
    """ % variant_id
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0])) for r in result])


@app.route('/heterozygous_individuals_json/<variant_id>')
@requires_auth
def get_heterozygous_individuals(variant_id):
    s="""
    MATCH
    (v)-[:HetVariantToPerson]-(p:Person)
    WHERE v.variantId='%s'
    RETURN p
    """ % variant_id
    db_session = neo4j_driver.session()
    result=db_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0])) for r in result])





