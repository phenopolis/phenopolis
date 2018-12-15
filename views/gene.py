from views import *
from lookups import *
from orm import *
import rest as annotation
import requests
from config import config
from vcf import vcf_query
import hashlib
from bson.json_util import dumps

'''
defs
'''
def hide_id_for_demo(data):
    if not data: return
    for k,v in data['patients'].items():
        # hide hpo
        v['hpo'] = ['hidden']
        # hide variants
        v['variants'] = ['hidden_'+hashlib.sha224(i).hexdigest()[:6] for i in v['variants']]
        # hide p_id
        new_p = 'hidden_'+hashlib.sha224(k).hexdigest()[:6]
        data['patients'][new_p] = data['patients'].pop(k)

    for k1,v1 in data['data'].items():
        for k2,v2 in v1['p'].items():
            v1['p'][k2] = ['hidden_'+hashlib.sha224(i).hexdigest()[:6] for i in v2]

    for k,v in data['variants'].items():
        new_v = 'hidden_'+hashlib.sha224(k).hexdigest()[:6]
        data['variants'][new_v] = data['variants'].pop(k)

'''
routes
'''
@app.route('/gene/<gene_id>',methods=['GET'])
#@cache.cached(timeout=24*3600)
@requires_auth
def gene_page(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'):
        gene=db.genes.find_one({'gene_name': gene_id}, projection={'_id': False})
        #if not gene: gene=db.genes.find_one({'other_names': gene_id}, projection={'_id': False})
        if not gene: return gene_id+' does not exist'
        gene_id=gene['gene_id']
    else:
        gene=db.genes.find_one({'gene_id':gene_id})
        if not gene: return gene_id+' does not exist'
    if session['user'] == 'demo' and gene_id not in ['ENSG00000156171','ENSG00000119685']: return 'Sorry you are not permitted to see these genes in demo account, please contact us to setup an account!'
    variants=db.variants.find({'genes':gene_id},projection={'_id':False})
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
    gene_hpo = db.gene_hpo.find_one({'gene_id':gene_id},{'_id':0})
    patients_status = {}
    if session['user'] == 'demo': hide_id_for_demo(gene_hpo) 
    else:
    # get patients status, solved? candidate genes? Only work when user is not demo for the time-being. Will probably change data struture later on to make it work for demo too
        all_patients = gene_hpo['patients'].keys()
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
    pli=get_db('exac').pli.find_one({'gene':gene['gene_name']})
    if pli:
        pli=pli['pLI']
    else:
        pli=-1
    return render_template('gene.html', 
            title=gene['gene_name'],
            gene=gene,
            pli=pli,
            table_headers=table_headers,
            phenogenon = json.dumps(gene_hpo) if gene_hpo else {},
            simreg = simreg,
            individuals=individuals,
            hpo_terms_json = json.dumps(hpo_terms),
            patients_status = dumps(patients_status),
            hpo_terms=hpo_terms_dict)


@app.route('/gene_json/<gene_id>',methods=['GET','POST'])
@requires_auth
def gene_json(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene=db.genes.find_one({'gene_id':gene_id})
    del gene['_id']
    variants=db.variants.find({'genes':gene_id})
    return json.dumps(gene)


# get sequence given region, and highlight the region. useful for design primers
@app.route('/sequence')
@requires_auth
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
    if not r.ok: return r.raise_for_status()
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
    return render_template('primer3_sequence.html', result = result,)
    
@app.route('/gene_phenogenon_json/<gene_id>',methods=['GET','POST'])
@requires_auth
def gene_phenogenon_json(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene=db.gene_hpo_new.find_one({'gene_id':gene_id})
    del gene['_id']
    return json.dumps(gene)

@app.route('/gene/',methods=['GET'])
@requires_auth
def gene():
   gene_id=request.args.get('id')
   x=json.loads(file('tests/data/TTLL5.json','r').read())
   return json.dumps(x)
    
# test
@app.route('/test')
def test():
    '''
    random test
    '''
    retnet_f=app.config['RETNET_JSON']
    RETNET = json.load(open(retnet_f, 'r'))
    relations = []
    genes = []
    omims = []
    for g, v in RETNET.iteritems():
        genes.append(g)
        for o in v['omim']:
            omims.append(o)
            relations.extend([(g,o)])
    return render_template('test.html', relations = json.dumps(relations), genes = json.dumps(list(set(genes))), omims = json.dumps(list(set(omims))))


    
