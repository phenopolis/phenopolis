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
    if not lookup_patient(db=get_db(app.config['DB_NAME_USERS']),user=session['user'],external_id=individual): return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual})
    if patient2 is None or 'report_id' not in patient2 or patient is None: return 'patient not loaded'
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
    return jsonify(result)



@app.route('/individual2/<individual>')
@requires_auth
def individual_page2(individual):
    # make sure that individual is accessible by user
    if not lookup_patient(db=get_db(app.config['DB_NAME_USERS']),user=session['user'],external_id=individual): return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual})
    if patient2 is None:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
        return redirect(referrer+'/load_individual/'+individual)
    patient['report_id']=patient2['report_id']
    patient['features']=patient2.get('features',[])
    patient['sex'] = patient2['sex']
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
    # get pubmedbatch table
    pubmedbatch = patient.get('pubmedbatch',{})
    # candidate genes
    patient['genes'] = patient2.get('genes',[])
    # solved genes
    patient['solved'] = patient2.get('solved',[])
    genes = {}
    # is this still updating?
    update_status = pubmedbatch.get('status', 0);
    # get known and retnet genes
    known_genes = open(app.config['RETNET_KNOWN_GENES'], 'r').readline().strip().split()
    RETNET  = json.load(open(app.config['RETNET_JSON'], 'r'))
    # sort the table, and add known / retnet gene annotation
    #for var_type in ['rare_homozygous', 'compound_hets', 'rare_variants']:
    for var_type in []:
        # mark the genes for retnet
        for ge in pubmedbatch[var_type]['result']:
            if ge['HUGO'] in known_genes:
                ge['ref(pubmedID)'][0]['known'] = 1
                # give the rest a minimal score to keep them on the list
                ge['ref(pubmedID)'][0]['total_score'] = max(1, ge['ref(pubmedID)'][0]['total_score'])
            else:
                ge['ref(pubmedID)'][0]['known'] = 0
            if ge['HUGO'] in RETNET:
                ge['ref(pubmedID)'][0]['disease'] = RETNET[ge['HUGO']]['disease']
                ge['ref(pubmedID)'][0]['omim'] = RETNET[ge['HUGO']]['omim']
                ge['ref(pubmedID)'][0]['mode'] = RETNET[ge['HUGO']]['mode']
                # reassign total score according to mode
                if RETNET[ge['HUGO']]['mode'] == 'd' or RETNET[ge['HUGO']]['mode'] == 'x':
                    ge['ref(pubmedID)'][1] = max(100, ge['ref(pubmedID)'][1])
                elif var_type ==  'rare_variants':
                    # not searching doimant, also assign others to 100
                    ge['ref(pubmedID)'][1] = max(100, ge['ref(pubmedID)'][1])
                else:
                    # give the rest a minimal score to keep them on the list
                    ge['ref(pubmedID)'][1] = max(1, ge['ref(pubmedID)'][1])
                # add pubmed result
                #this['ref(pubmedID)'] = [genes[gene_name], genes[gene_name]['total_score']]
        # now sort the table on the new scores        
        pubmedbatch[var_type]['result'] = sorted(pubmedbatch[var_type]['result'], key=lambda k: k['ref(pubmedID)'][1], reverse=True)
        # get genes for different var_type
        #variants = [p['variant_id'] for p in patient[var_type]]
        genes[var_type] = [g['HUGO'] for g in pubmedbatch[var_type]['result']]
        # delete pubmed_score column, if there is one. 
        if 'pubmed_score' in pubmedbatch[var_type]['header']:
            del pubmedbatch[var_type]['header'][pubmedbatch[var_type]['header'].index('pubmed_score')]
    # get combinatorics of features to draw venn diagram
    feature_combo = []
    feature_venn = []
    for i in range(len(hpo_terms)):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    #venn_ind = -1
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
    print feature_venn
    gene_info=dict()
    for v in patient['rare_variants']:
        if 'HUGO' not in v: v['HUGO']=''
        gene=v['HUGO'].upper() 
        gene_info[gene]=dict()
        if gene in known_genes: gene_info[gene]['known']=True
        if gene not in RETNET: continue
        gene_info[gene]['disease'] = RETNET[gene]['disease']
        gene_info[gene]['omim'] = RETNET[gene]['omim']
        gene_info[gene]['mode'] = RETNET[gene]['mode']
    genes['homozygous_variants']=[v.get('HUGO','').upper() for v in patient['homozygous_variants']]
    genes['compound_hets']=[v.get('HUGO','').upper() for v in patient['compound_hets']]
    genes['rare_variants']=[v.get('HUGO','').upper() for v in patient['rare_variants']]
    genes_pubmed=dict()
    for v in patient['rare_variants']:
        hugo=v['HUGO']
        genes_pubmed[hugo]=get_db('pubmedbatch').cache.find_one( {'key':re.compile(hugo+'[_ ].*')} )
    # figure out the order of columns from the variant row
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/variant_row.tmpl','r').read())
    print table_headers
    return render_template('individual.html', 
            external_id = individual,
            patient=patient,
            table_headers=table_headers,
            pubmedbatch=pubmedbatch,
            pubmed_db=get_db('pubmed_cache'),
            features = hpo_terms,
            genes = genes,
            hpo_gene = hpo_gene,
            gene_info=gene_info,
            genes_pubmed = genes_pubmed,
            update_status = update_status,
            feature_venn = feature_venn)


