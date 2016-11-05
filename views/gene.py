from views import *
from lookups import *
from orm import *
import rest as annotation
import requests
import primer3
import myvariant
from vcf import vcf_query
import hashlib
from bson.json_util import dumps
'''
defs
'''
def hide_hpo_for_demo(data):
    for mode in ['hom_comp','het']:
        for k1,v1 in data[mode].iteritems():
            # hide hpo
            for k2,v2 in v1['data'].iteritems():
                v2['hpo'] = ['hidden']
            # hide p_id
            for k2 in v1['data'].keys():
                v1['data']['hidden_'+hashlib.sha224(k2).hexdigest()[:6]] = v1['data'].pop(k2)
'''
routes
'''
@app.route('/gene/<gene_id>',methods=['GET'])
@requires_auth
def gene_page(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene=db.genes.find_one({'gene_id':gene_id})
    variants=db.variants.find({'genes':gene_id})
    gene['variants']=[Variant(variant_id=v['variant_id'],db=db) for v in variants]
    individuals=dict()
    for v in gene['variants']:
        v.canonical_hgvs=dict(zip( v.canonical_hgvsp, v.canonical_hgvsc))
        v.__dict__['protein_mutations']=dict([(p,p.split(':')[1],) for p in v.canonical_hgvsp if ':' in p])
        for s in v.het_samples:
            if v.HET_COUNT < 10:
                individuals[s]=individuals.get(s,[])+[v]
    print(gene['gene_id'])
    hpo_terms=hpo_db.gene_hpo.find_one({'gene_id':gene['gene_id']})
    if hpo_terms:
        #hpo_terms=dict(zip(hpo_terms['hpo_terms'],hpo_terms['hpo_names']))
        hpo_terms=hpo_terms['hpo_terms']
    else:
        hpo_terms=hpo_db.genes_pheno.find_one({'gene':gene['gene_name']})
        if hpo_terms:
            hpo_terms=hpo_terms['hpo']
        else:
            hpo_terms=[]
    hpo_terms_dict=dict()
    for hpo_id in hpo_terms:
        hpo_terms_dict[hpo_id]=hpo_db.hpo.find_one({'id':hpo_id})
    gene_hpo = db.gene_hpo.find_one({'gene_id':gene_id})
    patients_status = {}
    if session['user'] == 'demo': hide_hpo_for_demo(gene_hpo) 
    else:
    # get patients status, solved? candidate genes? Only work when user is not demo for the time-being. Will probably change data struture later on to make it work for demo too
        all_patients = frozenset(gene_hpo['het'].get('HP:0000001',{'data':{}})['data'].keys()) | frozenset(gene_hpo['hom_comp'].get('HP:0000001',{'data':{}})['data'].keys())
        patients_status = dict([(i['external_id'],i) for i in patient_db.patients.find({'external_id':{'$in':list(all_patients)}},{'external_id':1,'solved':1,'genes':1})])
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/gene-page-tabs/gene_variant_row.tmpl','r').read())
    # get simreg
    simreg_data = list(db.simreg.find({'gene':gene_id}))
    simreg = {'rec':{'data':[],'p':None},'dom':{'data':[],'p':None}}
    for mode in ['rec','dom']:
        temp = [i for i in simreg_data if i['mode'] == mode]
        if not temp: continue
        simreg[mode]['p'] = temp[0]['p']
        # convert it to array
        simreg[mode]['data'] = temp[0]['phi'].values()
        # sort desc
        simreg[mode]['data'] = sorted(simreg[mode]['data'], key=lambda x: x['prob'], reverse=True)
    return render_template('gene.html', 
            gene=gene,
            table_headers=table_headers,
            dot_hom_comp = json.dumps(gene_hpo['hom_comp']),
            dot_het = json.dumps(gene_hpo['het']),
            simreg = simreg,
            individuals=individuals,
            hpo_terms_json = json.dumps(hpo_terms),
            patients_status = dumps(patients_status),
            hpo_terms=hpo_terms_dict)

@app.route('/gene2/<gene_id>',methods=['GET'])
#@requires_auth
def gene_page2(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    hpo=request.args.get('hpo')
    hpo_freq = get_hpo_size_freq('hpo_freq_2016-7.tsv')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    hom_comp_file = os.path.join('dot',gene_id+'_hom_comp.json')
    het_file = os.path.join('dot',gene_id+'_het.json')
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    # get lof force
    # lof_p_hpo = get_lof_p_hpo(gene_id, db, patient_db)
    # rare_p_hpo = get_rare_var_p_hpo(gene_id, db, patient_db)
    # force_test = json.dumps(draw_force_graph(rare_p_hpo['het'], hpo_freq, hpo_db))
    hom_inf = open(hom_comp_file, 'r')
    het_inf = open(het_file,'r')
    hom_comp = json.load(hom_inf)
    het = json.load(het_inf)
    print '======'
    print session['user']
    print '====='
    #dot_hom_comp = json.dumps(hide_hpo_for_demo(hom_comp)) if session['user'] == 'demo' else json.dumps(hom_comp)
    dot_hom_comp = json.dumps(hom_comp)
    #dot_het = json.dumps(hide_hpo_for_demo(het)) if session['user'] == 'demo' else json.dumps(het)
    dot_het = json.dumps(het)
    return render_template('gene2.html', 
            gene_id = gene_id,
            gene_name = gene_name,
            dot_hom_comp = dot_hom_comp,
            dot_het = dot_het,
            hpo_freq = json.dumps(hpo_freq))

# get sequence given region, and highlight the region. useful for design primers
@app.route('/sequence')
def sequence():
    var_id = request.args.get('variant_id')
    symbol = request.args.get('symbol')
    paddings = int(request.args.get('paddings') or 100)
    build = request.args.get('build') or 'grch37'
    (chrom,start,ref) = var_id.split('-',3)[:3]
    start = int(start)
    end = max(start + len(ref) - 1, start)
    strand = request.args.get('strand') or '1'
    margins = int(request.args.get('margins') or 500)
    gc_opt = float(request.args.get('gc_opt') or 50)
    gc_min = float(request.args.get('gc_min') or 35)
    gc_max = float(request.args.get('gc_max') or 65)
    primer_size_opt = int(request.args.get('primer_size_opt') or 20)
    primer_size_min = int(request.args.get('primer_size_min') or 18)
    primer_size_max = int(request.args.get('primer_size_max') or 25)
    primer_tm_opt = float(request.args.get('primer_tm_opt') or 59)
    primer_tm_min = float(request.args.get('primer_tm_min') or 55)
    primer_tm_max = float(request.args.get('primer_tm_max') or 67)
    PCR_size_range = request.args.get('PCR_size_range') or '100-400'
    SEQUENCE_EXCLUDED_REGION = []

    # get sequence
    server = "http://%s.rest.ensembl.org" % build
    ext = '''/sequence/region/human/%(chrom)s:%(start)s..%(end)s:%(strand)s?expand_5prime=%(margins)s;expand_3prime=%(margins)s;''' % locals()
    r = requests.get(server+ext, headers={'Content-Type':'application/json' })
    if not r.ok:
        return r.raise_for_status()
    decoded = r.json()

    # run primer3 to get primers suggestion
    # useful keys:
    # PRIMER_LEFT_0: (start(0 based), length)
    # PRIMER_LEFT_0_SEQUENCE
    # PRIMER_LEFT_0_TM
    # PRIMER_LEFT_0_GC_PERCENT
    # PRIMER_PAIR_0_PRODUCT_SIZE

    # region sets the paddings to include in the sequencing.Default with 100 bp on each side.
    region = [margins - paddings, end - start + 2*paddings ]
    seq = str(decoded['seq'])
    primer =  primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': 'phenopolis',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': region
        },
        {
            'PRIMER_OPT_SIZE': primer_size_opt,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': primer_size_min,
            'PRIMER_MAX_SIZE': primer_size_max,
            'PRIMER_OPT_TM': primer_tm_opt,
            'PRIMER_MIN_TM': primer_tm_min,
            'PRIMER_MAX_TM': primer_tm_max,
            'PRIMER_MIN_GC': gc_min,
            'PRIMER_OPT_GC': gc_opt,
            'PRIMER_MAX_GC': gc_max,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE':[int(PCR_size_range.split('-')[0]),
                                         int(PCR_size_range.split('-')[1]) ]
        })
    if 'PRIMER_RIGHT_0' not in primer:
        # return 'Cannot pick any primers for the given sequence'
        left=left_tm=left_gc=right=right_tm=right_gc=product_length='NA'
        seq = seq[:margins-1] + '<span class="highlight">' + seq[margins-1:margins+end-start] + '</span>' + seq[margins+end-start:]
    else:
    # formulate sequence
        seq = seq[:primer['PRIMER_RIGHT_0'][0]-primer['PRIMER_RIGHT_0'][1]+1] + '<span class="primer">' + seq[primer['PRIMER_RIGHT_0'][0]-primer['PRIMER_RIGHT_0'][1]+1:primer['PRIMER_RIGHT_0'][0]+1] + '</span>' + seq[primer['PRIMER_RIGHT_0'][0]+1:]
        seq = seq[:margins] + '<span class="highlight">' + seq[margins:margins+end-start+1] + '</span>' + seq[margins+end-start+1:]
        seq = seq[:primer['PRIMER_LEFT_0'][0]] + '<span class="primer">' + seq[primer['PRIMER_LEFT_0'][0]:primer['PRIMER_LEFT_0'][0]+primer['PRIMER_LEFT_0'][1]] + '</span>' + seq[primer['PRIMER_LEFT_0'][0]+primer['PRIMER_LEFT_0'][1]:]
        left=primer['PRIMER_LEFT_0_SEQUENCE']
        left_tm=primer['PRIMER_LEFT_0_TM']
        left_gc=primer['PRIMER_LEFT_0_GC_PERCENT']
        right=primer['PRIMER_RIGHT_0_SEQUENCE']
        right_tm=primer['PRIMER_RIGHT_0_TM']
        right_gc=primer['PRIMER_RIGHT_0_GC_PERCENT']
        product_length=primer['PRIMER_PAIR_0_PRODUCT_SIZE']
    
    # construct object
    result = {'seq': seq,
            'var_id': var_id,
            'symbol': symbol,
            'left': left,
            'left_tm': left_tm,
            'left_gc': left_gc, 
            'right': right,
            'right_tm': right_tm,
            'right_gc': right_gc,
            'product_length': product_length,
            'paddings': paddings,
            'margins': margins,
            'gc_opt': gc_opt,
            'gc_max': gc_max,
            'gc_min': gc_min,
            'primer_tm_opt': primer_tm_opt,
            'primer_tm_max': primer_tm_max,
            'primer_tm_min': primer_tm_min,
            'primer_size_opt': primer_size_opt,
            'primer_size_min': primer_size_min,
            'primer_size_max': primer_size_max,
            'PCR_size_range': PCR_size_range,
            }

    return render_template('primer3_sequence.html',
                            result = result,
                            )
    

# test
@app.route('/test')
def test():
    '''
    random test
    '''
    retnet_f = 'retnet.json'
    RETNET = json.load(open(retnet_f, 'r'))
    relations = []
    genes = []
    omims = []
    for g, v in RETNET.iteritems():
        genes.append(g)
        for o in v['omim']:
            omims.append(o)
            relations.extend([(g,o)])

    return render_template('test.html',
                           relations = json.dumps(relations),
                           genes = json.dumps(list(set(genes))),
                           omims = json.dumps(list(set(omims)))
                           )


@app.route('/gene_json/<gene_id>',methods=['GET','POST'])
def gene_json(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    print(gene_name)
    hpo_string=lookups.get_gene_hpo(get_db('hpo'),gene_name)
    print(hpo_string)
    variants=db.variants.find({'genes': gene_id}, fields={'_id': False})
    #if gene_id in app.config['GENES_TO_CACHE']:
        #return open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id))).read()
    #else:



