from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import orm


@app.route('/individual_json/<individual>')
@requires_auth
def individual_json(individual):
    db=get_db()
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual},{'_id':False})
    print(patient2)
    return jsonify(result=patient2)

@app.route('/individual/<individual>')
@requires_auth
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
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual})
    print patient2
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
    known_genes = open('ret_known_genes.txt', 'r').readline().strip().split()
    RETNET  = json.load(open('retnet.json', 'r'))
    # get combinatorics of features to draw venn diagram
    feature_combo = []
    feature_venn = []
    print len(hpo_terms)
    for i in range(len(hpo_terms[1:5])):
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
    gene_info=dict()
    individuals=dict()
    for v in patient['rare_variants']:
        gene=v['canonical_gene_name_upper']
        gene_info[gene]=dict()
        if gene in known_genes: gene_info[gene]['known']=True
        if gene not in RETNET: continue
        gene_info[gene]['disease'] = RETNET[gene]['disease']
        gene_info[gene]['omim'] = RETNET[gene]['omim']
        gene_info[gene]['mode'] = RETNET[gene]['mode']
        if 'het_samples' not in v: print(v)
        for s in v['het_samples']:
            if v['HET_COUNT'] < 10:
                individuals[s]=individuals.get(s,[])+[v]
    genes['homozygous_variants']=[v['canonical_gene_name_upper'] for v in patient['homozygous_variants']]
    genes['compound_hets']=[v['canonical_gene_name_upper'] for v in patient['compound_hets']]
    genes['rare_variants']=[v['canonical_gene_name_upper'] for v in patient['rare_variants']]
    genes_pubmed=dict()
    for v in patient['rare_variants']+patient['homozygous_variants']+patient['compound_hets']:
        g=v['canonical_gene_name_upper']
        # gene_id is used to get gene-hpo analysis result
        temp = lookups.get_gene_by_name(get_db(), g)
        v['gene_id'] = temp['gene_id'] if temp else None
        genes_pubmed[g]=get_db('pubmedbatch').cache.find_one( {'key':g+'_'+patient.get('pubmed_key','')} )
        # print(g, genes_pubmed[g])
    # figure out the order of columns from the variant row
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/individual-page-tabs/individual_variant_row.tmpl','r').read())
    print table_headers
    # get a list of genes related to retinal dystrophy. only relevant to subset group of ppl. talk to Jing or Niko for other cohorts. Note that dominant p value only counts paitents with 1 qualified variant on the gene. 
    # current setting: unrelated, exac_af 0.01 for recessive, 0.001 for dominant, cadd_phred 15
    retinal_genes = {}
    if individual[:4] == 'IRDC':
        # retinal dystrophy == HP:0000556
        retinal_genes_raw = db.hpo_gene.find_one({'hpo_id':'HP:0000556'})
        # transform the data for easy access
        for mode in ['recessive','dominant']:
            retinal_genes[mode] = dict([(i['gene_id'],i['p_val']) for i in retinal_genes_raw['data']['unrelated'][mode]])
    return render_template('individual.html', 
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
            genes_pubmed = genes_pubmed,
            update_status = update_status,
            retinal_genes = retinal_genes,
            feature_venn = feature_venn)



@app.route('/individual_update/<individual>')
def individual_update(individual):
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    print 'UPDATE'
    print p
    print get_db('patients').patients.update({'external_id':individual},{'$set':p},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)



