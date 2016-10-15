from views import *
from load_individual import load_patient 

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt


#def get_db(dbname=None):
#    """
#    Opens a new database connection if there is none yet for the
#    current application context.
#    """
#    if dbname is None: dbname=app.config['DB_NAME']
#    if not hasattr(g, 'db_conn'):
#        g.db_conn=dict()
#        g.db_conn[dbname] = connect_db(dbname)
#    elif dbname not in g.db_conn:
#        g.db_conn[dbname] = connect_db(dbname)
#    return g.db_conn[dbname]



@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'): g.autocomplete_strings = [s.strip() for s in open(os.path.join(app.config['UCLEX_FILES_DIRECTORY'], 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = str(request.args.get('query'))
    #for n in dir(request): print(n, getattr(request,n))
    #print(request.HTTP_REFERER)
    print(request.referrer)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    else:
        referrer=''
    #u.netloc
    print(referrer)
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('{}/gene/{}'.format(referrer,identifier))
    elif datatype == 'transcript':
        return redirect('{}/transcript/{}'.format(referrer,identifier))
    elif datatype == 'variant':
        return redirect('{}/variant/{}'.format(referrer,identifier))
    elif datatype == 'region':
        return redirect('{}/region/{}'.format(referrer,identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('{}/dbsnp/{}'.format(referrer,identifier))
    elif datatype == 'hpo':
        return redirect('{}/hpo/{}'.format(referrer,identifier))
    elif datatype == 'mim':
        return redirect('{}/mim/{}'.format(referrer,identifier))
    elif datatype == 'error':
        return redirect('{}/error/{}'.format(referrer,identifier))
    elif datatype == 'not_found':
        return redirect('{}/not_found/{}'.format(referrer,identifier))
    else:
        raise Exception


@app.route('/patient/<patient_str>')
def get_patient(patient_str):
    pass

@app.route('/variant3/<variant_str>')
def variant_page3(variant_str):
    db=get_db()
    variant=db.variants.find_one({'VARIANT_ID':variant_str})
    patients=[p for p in db.patients.find({'external_id':{'$in': variant['HET']+variant['HOM']}})]
    hpo_terms=[p['features'] for p in patients]
    print(hpo_terms)
    print 'Rendering variant: %s' % variant_str
    return render_template( 'test.html', variant=variant)



@app.route('/variant/<variant_str>')
def variant_page(variant_str):
    db = get_db()
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    pos = int(pos)
    # pos, ref, alt = get_minimal_representation(pos, ref, alt)
    xpos = get_xpos(chrom, pos)
    variant = lookups.get_variant(db, xpos, ref, alt)
    print(variant)
    if variant is None:
        variant = {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt
        }
    consequences = None
    ordered_csqs = None
    if 'vep_annotations' in variant:
        variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
        ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
        ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
        consequences = defaultdict(lambda: defaultdict(list))
        for annotation in variant['vep_annotations']:
            annotation['HGVS'] = get_proper_hgvs(annotation)
            consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
    base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    any_covered = any([x['has_coverage'] for x in base_coverage])
    metrics = lookups.get_metrics(db, variant)
    # check the appropriate sqlite db to get the *expected* number of
    # available bams and *actual* number of available bams for this variant
    sqlite_db_path = os.path.join( app.config["READ_VIZ_DIR"], "combined_bams", chrom, "combined_chr%s_%03d.db" % (chrom, pos % 1000))
    print(sqlite_db_path)
    try:
        read_viz_db = sqlite3.connect(sqlite_db_path)
        n_het = read_viz_db.execute("select n_expected_samples, n_available_samples from t " "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'het')).fetchone()
        n_hom = read_viz_db.execute("select n_expected_samples, n_available_samples from t " "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'hom')).fetchone()
        read_viz_db.close()
    except Exception, e:
        logging.debug("Error when accessing sqlite db: %s - %s", sqlite_db_path, e)
        n_het = n_hom = None
    read_viz_dict = {
        'het': {'n_expected': n_het[0] if n_het is not None and n_het[0] is not None else -1, 'n_available': n_het[1] if n_het and n_het[1] else 0,},
        'hom': {'n_expected': n_hom[0] if n_hom is not None and n_hom[0] is not None else -1, 'n_available': n_hom[1] if n_hom and n_hom[1] else 0,},
    }
    for het_or_hom in ('het', 'hom',):
        #read_viz_dict[het_or_hom]['some_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] > 0)    and (read_viz_dict[het_or_hom]['n_expected'] - read_viz_dict[het_or_hom]['n_available'] > 0)
        read_viz_dict[het_or_hom]['all_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] != 0) and (read_viz_dict[het_or_hom]['n_available'] == 0)
        read_viz_dict[het_or_hom]['readgroups'] = [ '%(chrom)s-%(pos)s-%(ref)s-%(alt)s_%(het_or_hom)s%(i)s' % locals() for i in range(read_viz_dict[het_or_hom]['n_available']) ]   #eg. '1-157768000-G-C_hom1', 
        read_viz_dict[het_or_hom]['urls'] = [ os.path.join('combined_bams', chrom, 'combined_chr%s_%03d.bam' % (chrom, pos % 1000)) for i in range(read_viz_dict[het_or_hom]['n_available']) ]
    print 'Rendering variant: %s' % variant_str
    return render_template(
        'variant.html',
        variant=variant,
        base_coverage=base_coverage,
        consequences=consequences,
        any_covered=any_covered,
        ordered_csqs=ordered_csqs,
        metrics=metrics,
        read_viz=read_viz_dict,
    )

# AJAX
# Not finished
@app.route('/chisqu/<variant_str>',methods=['GET','POST'])
def chisq(variant_str):
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    print(region)
    records=tb.fetch(region=region)
    geno=dict(zip(headers, [r.split('\t') for r in records][0]))
    samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
    #d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
    #headers=file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().replace('#','').split('\t')
    #d=csv.DictReader(file('/data/UCLpheno/uclex-hpo.txt','r'),delimiter='\t')
    res=jsonify(result=hpo_patients)
    return res


def stream_template(template_name, **context):
    app.update_template_context(context)
    t = app.jinja_env.get_template(template_name)
    rv = t.stream(context)
    rv.enable_buffering(5)
    return rv

@app.route('/my-large-page.html')
def render_large_template():
    rows = iter_all_rows()
    return Response(stream_template('the_template.html', rows=rows))



@app.route('/stream')
def streamed_response():
    def generate():
         yield 'Hello '
         yield request.args['name']
         yield '!'
    return Response(stream_with_context(generate()))

def generate_patient_table():
    def get_variants(variant_ids):
        for vid in db.variants.find({'variant_id':{'$in':variant_ids}}):
            yield 

'''
serve the Vincent annotated csv files
'''
@app.route('/download/send_csv', methods=['GET','POST'])
@requires_auth
def download_csv():
    conn=PhenotipsClient()
    p_id = request.args.get('p_id')
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=p_id,auth=auth)
    if not p: return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    folder = request.args.get('folder')
    path = '/slms/UGI/vm_exports/vyp/phenotips/DROPBOX/'+session['user']
    csv_file = os.path.join(path,folder, p_id + '.csv')
    filename = folder+'_'+p_id+'.csv'
    if not os.path.isfile(csv_file):
        return 'Oops, file not found!'
    return send_file(csv_file,
                     mimetype='text/csv',
                     attachment_filename=filename,
                     as_attachment=True)

@app.route('/individual/<individual>')
@requires_auth
def individual_page(individual):
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
    known_genes = open('ret_known_genes.txt', 'r').readline().strip().split()
    RETNET  = json.load(open('retnet.json', 'r'))
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

# shows each individual, 
# all_individuals
@app.route('/individuals')
@requires_auth
def individuals_page(page=None):
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
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'homozygous_variants_count':1,'compound_hets_count':1, 'rare_variants_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    patients=conn.get_patient(auth=auth).get('patientSummaries',[])
    eids=[p['eid'] for p in patients]
    print(eids)
    patients=get_db('patients').patients.find({'external_id':{'$in':eids}})
    #patients=get_db('patients').patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    individuals=[f(p) for p in patients if 'external_id' in p]
    # family_history":{"consanguinity":true}
    return render_template('individuals_page.html',individuals=individuals)


@app.route('/research_pubmed', methods=['POST'])
def research_pubmed():
    # use new search terms to update the individual-pubmedbatch table
    patient_id = request.form['p_id']
    search_term = request.form['OR']
    # update patient pubmed result status as running (1)
    db=get_db()
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch.status': 1}})
    # do the actual update
    #exit_status = subprocess.call(['python','individual_pubmedBatch.py', '-p', patient_id, '--OR', search_term])
    exit_status=0
    # reset update status to 0
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch.status': 0}})
    return str(exit_status)


@app.route('/hpo/<hpo_id>')
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
    print(len(hpo_ancestors))
    print([h['name'] for h in hpo_ancestors])
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
    pmids=[]
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
    return render_template('phenotype.html',hpo_id=hpo_id,hpo_name=hpo_name,individuals=[str(p['external_id']) for p in patients],genes=genes,pmids=pmids,variants=[])

# AJAX
# fetch patients iwth hpo term
@app.route('/fetch_hpo',methods=['GET','POST'])
def fetch_hpo():
    if request.method=='POST':
        hpo_ids=request.form['hpo_ids'].strip().split(',')
    else:
        hpo_ids=request.args.get('hpo_ids').strip().split(',')
    hpo_id=hpo_ids[0]
    print('HPO',hpo_id)
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    hpo_patients=[p['external_id'] for p in lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)]
    print('num patients',len(hpo_patients))
    res=jsonify(result=hpo_patients)
    return res

# AJAX
# fetch variants private to patients
# That is variants which are only seen in these patients and no one else.
@app.route('/fetch_private_variants',methods=['GET','POST'])
def fetch_private_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    db=get_db()
    if len(hpo_patients)==1:
        variants=db.variants.find({'PRIVATE_MUT':hpo_patients})
    else:
        #rsession=get_R_session()
        variants=rsession.r.private_variants(hpo_patients)
        #variants=[]
        print('private variants', variants)
        if type(variants) is str:
            variants=[variants]
        else:
            variants=variants.tolist()
    print('num of private variants',len(variants),)
    res=jsonify(result=variants)
    return res

# AJAX
# fetch common variants to patients
# That is variants which are seen in all these patients.
@app.route('/fetch_common_variants',methods=['GET','POST'])
def fetch_common_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    #rsession=get_R_session()
    #variants=rsession.r.common_variants(hpo_patients)
    variants=[]
    print('common variants', variants)
    if type(variants) is str:
        variants=[variants]
    else:
        variants=variants.tolist()
    print('num of common variants',len(variants),)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches variant record from db
@app.route('/fetch_variant',methods=['GET','POST'])
def fetch_variant():
    if request.method=='POST':
        variants=request.form['variants'].strip().split(',')
    else:
        variants=request.args.get('variants').strip().split(',')
    db=get_db()
    req_len=len(variants)
    variant_ids=map(lambda x: x.replace('_','-'),variants)
    variants=[v for v in db.variants.find({'variant_id':{'$in':variant_ids}}, fields={'_id': False})]
    ans_len=len(variants)
    print(req_len==ans_len)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches information from db
@app.route('/variant_count',methods=['GET','POST'])
def variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    #rsession=get_R_session()
    #res=jsonify(result={'variant_count':rsession.eval('sum(as.logical(variants[["%s"]]))' % external_id) , 'external_id':external_id})
    #return res

# AJAX
# fetches information from db
@app.route('/private_variant_count',methods=['GET','POST'])
def private_variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    db=get_db('patients')
    p=db.patients.find_one({'external_id':external_id})
    if 'PRIVATE_MUT' not in p: private_variant_count=0
    else: private_variant_count=len(p['PRIVATE_MUT'])
    res=jsonify(result={'variant_count': private_variant_count, 'external_id':external_id})
    return res


@app.route('/mim/<mim_id>')
def mim_page(mim_id):
    db=get_db('patients')
    print(str(mim_id))
    patients=[p for p in db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    patient_ids=[p['external_id'] for p in patients]
    print(phizz.query_disease([hpo_id]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    return render_template('test.html')

@app.route('/patient/<patient_id>')
def patient_page(patient_id):
    db=get_db()
    patients=[p for p in db.patients.find({'external_id': str(patient_id)})]
    print(patients)
    return None

@app.route('/Exomiser/<path:path>')
@requires_auth
def exomiser_page(path):
    #is this user authorized to see this patient?
    return send_from_directory('Exomiser', path)

@app.route('/example/')
@requires_auth
def example():
    return send_from_directory('templates', 'temp-plot.html')

@app.route('/gene/<gene_id>',methods=['GET'])
def gene_page(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    print(gene_name)
    hpo_string=lookups.get_gene_hpo(get_db('hpo'),gene_name)
    #if gene_id in app.config['GENES_TO_CACHE']:
        #return open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id))).read()
    #else:
    return get_gene_page_content(gene_id,hpo,hpo_string)


@app.route('/transcript2/<transcript_id>')
def transcript_page2(transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)
        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
        if t is None:
            gene = lookups.get_gene(db, transcript['gene_id'])
            gene['transcripts'] = lookups.get_transcripts_in_gene(db, transcript['gene_id'])
            variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
            coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
            t = render_template(
                'transcript.html',
                transcript=transcript,
                transcript_json=json.dumps(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=json.dumps(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=json.dumps(coverage_stats),
                gene=gene,
                gene_json=json.dumps(gene),
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)



@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    transcript = lookups.get_transcript(db, transcript_id)
    cache_key = 't-transcript-{}'.format(transcript_id)
    t = cache.get(cache_key)
    print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
    if t: return t
    variants=[v for v in db.variants.find({'Transcript':str(transcript_id)})]
    genes=list(set([variants['Gene'] for v in variants]))
    print(genes)
    cache.set(cache_key, t)
    return t



@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}'.format(region_id)
        t = cache.get(cache_key)
        print 'Rendering %sregion: %s' % ('' if t is None else 'cached ', region_id)
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                    csq_order=csq_order,
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
        return render_template(
            'region.html',
            rsid=rsid,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            csq_order=csq_order,
        )
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    )


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')

@app.route('/samples')
def samples_page():
    samples=pandas.read_csv('/slms/UGI/vm_exports/vyp/phenotips/HPO/hpo.txt')
    return render_template('samples.html',samples=samples.to_html(escape=False))


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.route('/read_viz/<path:path>')
def read_viz_files(path):
    full_path = os.path.abspath(os.path.join(app.config["READ_VIZ_DIR"], path))
    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(app.config["READ_VIZ_DIR"]):
        return "Invalid path: %s" % path
    logging.info("path: " + full_path)
    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(app.config["READ_VIZ_DIR"], path)
    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        logging.error(error_msg)
        return error_msg
    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset
    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)
    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))
    logging.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv


@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response


### all the mongodb reading/writing code


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_dbsnp_file, load_base_coverage, load_gene_models, load_constraint_information]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))
    [p.join() for p in all_procs]
    print('Done! Loading MNPs...')
    load_mnps()
    print('Done! Creating cache...')
    #create_cache()
    print('Done!')


def load_base_coverage():
    """ """
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation, e:
            print(e)
            # handle error when coverage_generator is empty
            pass  
    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.ensure_index('xpos')
    procs = []
    coverage_files = app.config['BASE_COVERAGE_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def load_variants_file():
    def load_variants(sites_file, i, n, db):
        for f in sites_file:
            print(f)
            variants_generator = parse_tabix_file_subset([f], i, n, get_variants_from_sites_vcf)
            try:
                db.variants.insert(variants_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when variant_generator is empty
    db = get_db('exac')
    db.variants.drop()
    print("Dropped db.variants")
    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')
    db.variants.ensure_index('variant_id')
    #sites_vcfs = app.config['SITES_VCFS']
    sites_vcfs=['/slms/UGI/vm_exports/vyp/phenotips/ExAC/0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz']
    print(sites_vcfs)
    #if len(sites_vcfs) > 1: raise Exception("More than one sites vcf file found: %s" % sites_vcfs)
    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    #pdb.set_trace()
    for i in range(num_procs):
        p = Process(target=load_variants, args=(sites_vcfs, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)


def load_constraint_information():
    db = get_db()
    db.constraint.drop()
    print 'Dropped db.constraint.'
    start_time = time.time()
    with gzip.open(app.config['CONSTRAINT_FILE']) as constraint_file:
        for transcript in get_constraint_information(constraint_file):
            db.constraint.insert(transcript, w=0)
    db.constraint.ensure_index('transcript')
    print 'Done loading constraint info. Took %s seconds' % int(time.time() - start_time)


def load_mnps():
    db = get_db()
    start_time = time.time()
    db.variants.ensure_index('has_mnp')
    print 'Done indexing.'
    while db.variants.find_and_modify({'has_mnp' : True}, {'$unset': {'has_mnp': '', 'mnps': ''}}): pass
    print 'Deleted MNP data.'
    with gzip.open(app.config['MNP_FILE']) as mnp_file:
        for mnp in get_mnp_data(mnp_file):
            variant = lookups.get_raw_variant(db, mnp['xpos'], mnp['ref'], mnp['alt'], True)
            db.variants.find_and_modify({'_id': variant['_id']}, {'$set': {'has_mnp': True}, '$push': {'mnps': mnp}}, w=0)
    db.variants.ensure_index('has_mnp')
    print 'Done loading MNP info. Took %s seconds' % int(time.time() - start_time)


def load_gene_models():
    db = get_db()
    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print 'Dropped db.genes, db.transcripts, and db.exons.'
    start_time = time.time()
    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript
    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)
    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)
    print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)
    # grab genes from GTF
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)
    print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)
    # and now transcripts
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)
    # Building up gene definitions
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)
    return []


def load_dbsnp_file():
    db = get_db()
    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty
        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)
    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']
    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"): num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())
    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)
    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


"""
Get the most recent common ancestor between two sets of hpo terms.
"""
def mrc_hpo():
    hpo_graph=get_hpo_graph()
    db=get_db()
    for var in db.variants.find():
        hpo_anc=[]
        for eid in list(set(var['HET']+var['HOM'])):
            patient=db.patients.find_one({'external_id':eid})
            if not patient: continue
            if 'features' not in patient: continue
            for f in patient['features']:
                fid=f['id']
                if not fid.startswith('HP'): continue
                hpo_anc.append(set(hpo_graph.get_ancestors(fid)))
        if not hpo_anc: continue
        if 'SYMBOL' not in var: continue
        var['ALL_HPO']=list(set(set.union(*hpo_anc)))
        var['SHARED_HPO']=list(set.intersection(*hpo_anc))
        print(var['VARIANT_ID'],var['SYMBOL'],len(var['HET']+var['HOM']),var['SHARED_HPO'],var['ALL_HPO'])
        db.variants.update({'VARIANT_ID':var['VARIANT_ID']},var,upsert=True)


#progressbar
'''
{
    'random_p_id':{
        'total':456,
        'count':123,
        'status':['running','done']
    },
    ...
}
'''
PROGRESS_BAR = {}

'''
initiate a progress instance
arg: total length of genes
return: progress_id
'''

def init_progress_bar(id,length):
    # check id
    if id in PROGRESS_BAR:
        if PROGRESS_BAR[id]['status'] != 'done':
            return 'the id already exists in PROGRESS_BAR'

    # initialise progress_bar
    PROGRESS_BAR[id] = {
            'total': length,
            'count':0,
            'message': '',
            'status':'running'
    }
    return 0

'''
update progress
arg: {
    id: id, 
    message: message,
    step: 1
    }
default step 1
'''

def update_progress_bar(obj):
    # check if id in PROGRESS_BAR
    if not obj['id'] in PROGRESS_BAR:
        return 'ID does not exist in PROGRESS_BAR'

    # update progress
    if not 'step' in obj:
        obj['step'] = 1
    PROGRESS_BAR[obj['id']]['count'] += obj['step']

    PROGRESS_BAR[obj['id']]['message'] = obj['message']
    # done?
    if PROGRESS_BAR[obj['id']]['count'] == PROGRESS_BAR[obj['id']]['total']:
        PROGRESS_BAR[obj['id']]['status'] = 'done'

'''
kill a progress
'''

def kill_progress_bar(key):
    if key in PROGRESS_BAR:
        del PROGRESS_BAR[key]


'''
to check if an iterable is empty
'''
def peek(iterable):
    try:
        first = next(iterable)
    except RuntimeError:
        return None
    except StopIteration:
        return None
    return first, itertools.chain([first], iterable)


'''
find the freaking PID, Title or Abstract no matter what!
'''
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



"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def scrutinise(obj):
    print obj['smashed_all']
    if obj['lag']:
        obj['lag'] = obj['lag']/3600/24 # convert it to days
            # need to update
        search_results = Entrez.read(Entrez.esearch(db='pubmed', term=obj['smashed_all'], reldate=obj['lag'], datetype='pdat', usehistory='y'))
    else:
        # just search
        search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=obj['smashed_all'], usehistory='y'))
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'total_score':0}
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed",
                                   restart=0,
                                   retmax=50,
                                   retmode="xml",
                                   webenv=search_results['WebEnv'],
                                   query_key=search_results['QueryKey']
                                   )
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print('Received error from server %s' % err)
            else:
                print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    if peek(record):
        # got something. let's do some calculation
        for r in record:
            # calculate score
            score = 0
            pid = str(find_item(r, 'PMID'))
            abstract_list = find_item(r, 'AbstractText')
            # parse abstract
            abstract = ''
            if abstract_list:
                for a in abstract_list:
                    if hasattr(a, 'attributes') and 'Label' in a.attributes:
                        abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                        abstract = abstract + a + '<br/>'
                    else:
                        abstract = abstract + a
        
            title = find_item(r, 'ArticleTitle')
            if title:
                score = score + len(obj['reg'].findall(title))
            if abstract:
                score = score + len(obj['reg'].findall(abstract))
        
            # add result to genes[gene_name]
            if score:
                results['results'].append({
                    'id': pid,
                    'title': title,
                    'abstract': abstract,
                    'score': score
                })
                results['total_score'] = results['total_score'] + score
    results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results

def get_pred_score(obj):
    # for the batch_pubmed route.
    # calculate the pred score
    # [D/A].each = 10, [P].each = 5, [C].each = 6, [T/B/N].each = -1. If there is a splicing/insertion/deletion event, the score is set as 1000. Not given is set as 0
    # ref: https://github.com/plagnollab/DNASeq_pipeline/blob/master/GATK_v2/filtering.md
    pred = 0
    if ('Func' in obj and re.search('splic', obj['Func'])) or ('ExonicFunc' in obj and re.search(r'stop|frame|del|insert', obj['ExonicFunc'])):
        pred = 1000;
    else:
        for k in obj:
            if re.search('Pred', k):
                if obj[k] == 'D' or obj[k] == 'A':
                    pred = pred + 10
                elif obj[k] == 'P':
                    pred = pred + 5
                elif obj[k] == 'C':
                    pred = pred + 6
                elif obj[k] == 'T' or obj[k] == 'B' or obj[k] == 'N':
                    pred = pred - 1
                else:
                    pass
    return pred;


@app.route('/plot/<gene>')
def plot(gene):
    #db = get_db()
    #var=db.variants.find_one({'VARIANT_ID':'3_8775295_C_T'})
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    x=[i for i, r, in enumerate(d)]
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    y=[-math.log10(float(r['HCM.chisq.p'])) for r in d]
    print(x)
    print(y)
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    #layout = dict( yaxis = dict( type = 'log', tickvals = [ 1.5, 2.53, 5.99999 ]), xaxis = dict( ticktext = [ "green eggs", "& ham", "H2O", "Gorgonzola" ], tickvals = [ 0, 1, 2, 3, 4, 5 ]))
    labels=[r['VARIANT_ID'] for r in d]
    layout = Layout( xaxis = dict( ticktext=labels, tickvals=x ), title="p-value plot" )
    #Layout( title="p-value plot")
    plotly.offline.plot({
        "data": [
                Scatter(
                    x=x,
                    y=y
                    )
                ],
        "layout": layout
        }, filename='genes/%s-pvalues.html' % (gene,), auto_open=False)
    return send_from_directory('genes', '%s-pvalues.html' % gene,)

""" JINJA2 filer """
def highlight(text, list, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    for l in list:
        # remove (.*), escape +?.*
        l = re.sub(r'\(.*\)', '', l)
        l = re.sub(r'\+','\\+',l)
        l = re.sub(r'\?','\\?',l)
        l = re.sub(r'\.','\\.',l)
        l = re.sub(r'\*','\\*',l)
        l = re.sub(r'\[.*\]','',l)
        l = re.sub(r'\\', '\\\\',l)
        words = l.split(',')
        for w in words:
            # wrap w with brackets to be a catch group
            text = re.sub(r'(\b%s\b)' % w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight'] = highlight

def highlight2(text, kw, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    # remove (.*), escape +?.*
    for w in kw:
        # wrap w with brackets to be a catch group
        text = re.sub(r'(%s)'%w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight2'] = highlight2


@app.route('/load_individual/<individual>')
@requires_auth
def load_individual(individual):
    filename='/slms/UGI/vm_exports/vyp/phenotips/DROPBOX/rare_variants/%s.csv' % individual
    auth='%s:%s' % (session['user'],session['password2'],)
    p = Process(target=load_patient, args=(filename,auth))
    p.start()
    return 'Loading %s...' % individual


