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
from load_individual import load_patient


@app.route('/individual_json/<individual>')
@requires_auth
def individual_json(individual):
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    #if not patient: return jsonify(result=None)
    patient = patient_db.patients.find_one({'external_id':individual},{'_id':False})
    patient['features']=[f for f in patient['features'] if f['observed']=='yes']
    print(patient)
    return jsonify(result=patient)

@app.route('/individual/<individual>')
@requires_auth
@cache.cached(timeout=24*3600)
def individual_page(individual):
    #print 'full_path', request.full_path
    #print  'url_root', request.url_root
    #if session['user']=='demo': individual=decrypt(str(individual))
    # make sure that individual is accessible by user
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    if not p: return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual})
    if patient2 is None or 'report_id' not in patient2 or patient is None:
        print 'no patient in mongo'
        referrer=request.referrer
        if referrer:
            u = urlparse(referrer)
            referrer='%s://%s' % (u.scheme,u.hostname,)
            if u.port: referrer='%s:%s' % (referrer,u.port,)
        else:
            referrer=''
        url=referrer+'/load_individual/'+individual
        url=request.url_root.replace('http','https')+'/load_individual/'+individual
        print(url)
        return redirect(url)
    patient['report_id']=patient2['report_id']
    patient['sex'] = patient2['sex']
    patient['features'] = patient2['features']
    patient['family_history'] = patient2.get('family_history',[])
    hpo_ids=[f['id'] for f in patient['features'] if f['observed']=='yes']
    # TODO
    # mode of inheritance in hpo terms: HP:0000005
    #print lookups.get_hpo_children(hpo_db, 'HP:0000005')
    patient['global_mode_of_inheritance']=patient2.get('global_mode_of_inheritance',None)
    # minimise it
    hpo_ids = lookups.hpo_minimum_set(hpo_db, hpo_ids)
    hpo_terms = [(i, hpo_db.hpo.find_one({'id':i})['name'][0]) for i in hpo_ids]
    # this has missing HPO ids. see IRDC_batch2_OXF_3001 and #HP:0000593
    hpo_gene=dict()
    for hpo_id,hpo_term, in hpo_terms:
        hpo_gene[hpo_id] = []
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            #gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
            hpo_gene[hpo_id]=hpo_gene.get(hpo_id,[])+[gene_name]
    for k in hpo_gene: hpo_gene[k]=list(frozenset(list(hpo_gene[k])))
    print '========'
    print hpo_gene
    print '========'
    # get pubmedbatch scores
    pubmedbatch = {}
    if patient.get('pubmed_key',None):
        genes = [v.get('canonical_gene_name_upper',None) for v in patient['rare_variants']+patient['homozygous_variants']+patient['compound_hets']]
        pubmed_keys = ['_'.join([g,patient['pubmed_key']]) for g in set(genes)]
        pubmedbatch = list(get_db('pubmedbatch').cache.find({'key':{'$in':pubmed_keys}},{'key':1,'score':1,'_id':0}))
        if pubmedbatch:
            pubmedbatch = dict([(i['key'],i.get('score',None)) for i in pubmedbatch])
    # candidate genes
    patient['genes'] = patient2.get('genes',[])
    # solved genes
    patient['solved'] = patient2.get('solved',[])
    genes = {}
    # is this still updating?
    if type(pubmedbatch) is dict:
        update_status = pubmedbatch.get('status', 0)
    else:
        update_status=0
    # get known and retnet genes
    known_genes=[x['gene_name'] for x in db.retnet.find()]
    RETNET = dict([(i['gene_name'],i) for i in db.retnet.find({},projection={'_id':False})])
    # get combinatorics of features to draw venn diagram
    feature_combo = []
    feature_venn = []
    print len(hpo_terms)
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    #venn_ind = -1
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = '","'.join([hpo_terms[i][1] for i in combo])
        dic_key = '"' + dic_key + '"'
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':frozenset(hpo_gene.get(x,""))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = feature_venn[-1]['value'] & frozenset(hpo_gene[hpo_terms[combo[ind]][0]])
    print 'get pubmed score and RETNET'
    gene_info=dict()
    individuals=dict()
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_hets']:
        for v in patient[mm]:
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
    genes['homozygous_variants']=[v['canonical_gene_name_upper'] for v in patient['homozygous_variants']]
    genes['compound_hets']=[v['canonical_gene_name_upper'] for v in patient['compound_hets']]
    genes['rare_variants']=[v['canonical_gene_name_upper'] for v in patient['rare_variants']]
    print 'get annotation'
    for v in patient['rare_variants']+patient['homozygous_variants']+patient['compound_hets']:
        g=v['canonical_gene_name_upper']
        # gene_id is used to get gene-hpo analysis result
        temp = lookups.get_gene_by_name(get_db(), g)
        v['gene_id'] = temp['gene_id'] if temp else None
        v['canonical_hgvs']=dict(zip( v.get('canonical_hgvsp',''), v.get('canonical_hgvsc','')))
        v['protein_mutations']=dict([(p,p.split(':')[1],) for p in v.get('canonical_hgvsp','') if ':' in p])
        if 'FILTER' not in v: v['FILTER']=v['filter']
        if 'ID' not in v: v['ID']=''
        # print(g, genes_pubmed[g])
    # figure out the order of columns from the variant row
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?.*>",file('templates/individual-page-tabs/individual_variant_row.tmpl','r').read())
    if session['user']=='demo': table_headers=table_headers[:-1]
    print table_headers
    # get a list of genes related to retinal dystrophy. only relevant to subset group of ppl. talk to Jing or Niko for other cohorts. Note that dominant p value only counts paitents with 1 qualified variant on the gene. 
    # current setting: unrelated, exac_af 0.01 for recessive, 0.001 for dominant, cadd_phred 15
    print 'get phenogenon genes'
    retinal_genes = {}
    patient['pubmed_key']=patient.get('pubmed_key','')
    if individual[:4] == 'IRDC' or individual in ['WebsterURMD_Sample_GV4344','WebsterURMD_Sample_IC16489','WebsterURMD_Sample_SJ17898','WebsterURMD_Sample_SK13768']:
        # retinal dystrophy == HP:0000556
        retinal_genes_raw = db.hpo_gene.find_one({'hpo_id':'HP:0000556'})
        # transform the data for easy access
        for mode in ['recessive','dominant']:
            retinal_genes[mode] = dict([(i['gene_id'],-math.log10(i['p_val'])) for i in retinal_genes_raw['data']['unrelated'][mode]])
    return render_template('individual.html', 
            title=individual,
            external_id = individual,
            patient=patient,
            table_headers=table_headers,
            pubmedbatch=pubmedbatch,
            pubmed_db=get_db('pubmed_cache'),
            features = hpo_terms,
            genes = genes,
            individuals=individuals,
            hpo_gene = hpo_gene,
            gene_info=gene_info,
            update_status = update_status,
            retinal_genes = retinal_genes,
            feature_venn = feature_venn)



@app.route('/individual_update/<individual>')
@requires_auth
def individual_update(individual):
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth,cache=False)
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

@app.route('/load_individual/<individual>')
@requires_auth
def load_individual(individual):
    auth='%s:%s' % (session['user'],session['password2'],)
    p = Process(target=load_patient, args=(individual,auth))
    p.start()
    return 'Loading %s...' % individual




