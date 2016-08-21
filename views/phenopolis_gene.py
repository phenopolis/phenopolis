from views import *
from lookups import *
import annotation
import requests
import primer3
'''
some defs
'''
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def get_venn(p_hpo):
    # useless
    #p_hpo = {'p_id':{hpo:[(id,term),(id,term)], exac_af:[], uclex_af:[]}}
    all_p_ids = f7(p_hpo.keys())
    # set combo
    p_combo = []
    venn = []
    for i in range(len(all_p_ids)):
        p_combo.extend(itertools.combinations(range(len(all_p_ids)), i+1))
    return p_combo
    for combo in p_combo:
        # get intersections
        intersect = set(p_hpo[all_p_ids[combo[0]]]['hpo'])
        key = []
        for ind in combo:
            key.append(all_p_ids[ind])
            intersect = intersect & set(p_hpo[all_p_ids[ind]]['hpo'])
        # convert immutable list to mutable list, to be recognised by js
        value = [ list(i) for i in intersect]
        size = len(value)
        value = json.dumps(value)
        key = '","'.join(key)
        key = '"' + key + '"'
        # set a venn area
        if value:
            venn.append({'key':key, 'size': size, 'value':value})
    return venn

def draw_force_graph(patients, hpo_freq, hpo_db):
    '''
    patients = {p_id:{'hpo':[(HP:001, hel), ...]}}
    result = '{
      "nodes":[
          {"name":"Myriel","group":1, "hpo":[["HP:0000345","hell"],["HP..]] },
          {"name":"Napoleon","group":1, "hpo": ..},
          {"name":"Mlle.Baptistine","group":1, "hpo": ..},
          {"name":"Mme.Magloire","group":1, "hpo": ..}
      ],
      "links":[{"source":1,"target":0,"value":1},
      {"source":2,"target":0,"value":8},
      {"source":3,"target":0,"value":10},
      ]
  }
  set all to group one for the time being
    '''
    # convert patients to array to conserve order
    for p in patients:
        patients[p]['hpo'] = [list(i) for i in patients[p]['hpo']]
    patients = [{'name':key, 'group':1, 'hpo':patients[key]['hpo']} for key in patients]
    result = {'nodes':patients, 'links':[]}
    p_combo = itertools.combinations(range(len(patients)), 2)
    for pair in p_combo:
        p1 = [i[0] for i in patients[pair[0]]['hpo']]
        p2 = [i[0] for i in patients[pair[1]]['hpo']]
        value = patients_hpo_nearest_similarity(p1, p2, hpo_freq, hpo_db)
        if not value:
            continue
        result['links'].append({
            'source':pair[0],
            'target':pair[1],
            'value':value
            })
    return result
def patients_hpo_nearest_similarity(p1, p2, hpo_freq, hpo_db, cutoff=0.5):
    # calculate patients similarity given hpo terms
    # p1 = [HPid, HPid]
    '''
    hpo similarity is calculated by Lin's similarity
    find the most similar hpo terms between p1 and p2, and report the
        similarity score
    if nearest common ancestor's freq >= cutoff, no link
    '''
    # minimise
    p1 = hpo_minimum_set(hpo_db, hpo_ids=p1)
    p2 = hpo_minimum_set(hpo_db, hpo_ids=p2)
    result = 0.
    for i in p1:
        for j in p2:
            result = max( result, lin_similarity(i,j,hpo_freq,hpo_db, cutoff) )
    return result
    
def patients_hpo_average_similarity(p1, p2, hpo_freq, hpo_db, cutoff=0.5):
    # calculate patients similarity given hpo terms
    #p1=[HPid, HPid]
    '''
    hpo similarity is calculated by Lin's similarity
    sum the similarities and then normalised by the number of pairs
    if nearest common ancestor's freq >= cutoff, no link
    '''
    # minimise
    p1 = hpo_minimum_set(hpo_db, hpo_ids=p1)
    p2 = hpo_minimum_set(hpo_db, hpo_ids=p2)
    num_connections = len(p1) * len(p2)
    result = 0.
    for i in p1:
        for j in p2:
            result += lin_similarity(i,j,hpo_freq,hpo_db, cutoff)
    return result/num_connections

def IC(freq):
    # calculate information content
    return -math.log(freq)

def lin_similarity(h1, h2, hpo_freq, hpo_db, cutoff):
    '''
    hpo similarity is calculated by Lin's similarity
    s(t1,t2) = [ 2 * max(IC(t))] / (IC(t1) + IC(t2))
    where t is from common ancestors of t1 and t2: t = set(ancestor(t1)) & set(ancestor(t2))
    IC(t) = -log(frequency(t))
    '''
    t1 = hpo_freq[h1] 
    t2 = hpo_freq[h2]
    # nearest ancestor, and get freq
    n_a = get_hpo_nearest_common_ancestors(hpo_db, h1, h2, hpo_freq)
    if hpo_freq[n_a[0]] >= cutoff:
        return 0
    max_an_IC = IC(hpo_freq[n_a[0]])
    return max_an_IC * 2 / (IC(t1) + IC(t2))


def get_lof_p_hpo(gene_id, db, patient_db):
    # since af is not very important here, leave them as 0 for the time being.
    #return {'hom_comp':{p_id1: {hpo:[(HP:1234,hell)], exac_af:[0]},
    #        'het':{p_id2: [(HP:2345,yeah)]}}
    lof_gene = db.lof.find_one({'gene_id':gene_id})
    patients_hpo = {'hom_comp':{},'het':{}}
    if not lof_gene:
        return patients_hpo
    patients = lof_gene['patient_ids']
    for p in patients:
        hpo = get_patient_observed_hpo(p, patient_db)
        if len(patients[p]) > 1:
            patients_hpo['hom_comp'][p]= {'hpo': hpo, 'exac_af': [0], 'uclex_af': [0]}
        else:
            patients_hpo['het'][p]= {'hpo': hpo, 'exac_af': [0], 'uclex_af': [0]}
    return patients_hpo

def get_rare_var_p_hpo(gene_id, db, patient_db):
    #return {'hom_comp':{p_id1:{hpo: [(HP:1234,hell)], exac_af:[0.0021], uclex_af:[0.001]},
    #        'het':{p_id2: {hpo:[(HP:2345,yeah)], exac_af:[0.001,0.002],uclex_af:[0.002,0.001]}}

    # sometimes variant is not in vcf. move it to debug/bad_variants for inspection and later clean
    bad_var_file = open('views/debug/bad_variants', 'w')
    # get all variants on this gene
    all_vars = db.genes.find_one({'gene_id':gene_id})['variant_ids']
    results = {'hom_comp':{}, 'het':{}} 
    for v in all_vars:
        var = db.variants.find_one({'variant_id':v})
        exac_af = 0
        if var['in_exac']:
            if 'allele_freq' not in var['EXAC']:
                VAR = annotation.exac_anno(v)
                exac_af = VAR['allele_freq']
            else:
                exac_af = var['EXAC']['allele_freq']
        # not interested if af is > 0.01
        if float(exac_af) > 0.01:
            continue
        # get relevant info from vcf
        this = lookups.vcf_query(variant_str=v)
        if not this:
            bad_var_file.write(v+'\n')
            continue
        uclex_af = this['allele_freq']

        # dealing with hom patients. also add it to het. count hom as twice
        # will need to deal with both_het !!! their af are different!!!
        hom_p = this['hom_samples']
        for p in hom_p:
            populate_mode_p(results, 2, ['hom_comp', 'het'], p, exac_af, uclex_af, patient_db)

        # dealing with het patients. note to check length of exac_af. longer than one?
        # also added it to 'hom_comp'
        het = this['het_samples']
        for p in het:
            results['het'][p] = results['het'].get(p, {'exac_af':[], 'uclex_af':[]})
            modes = ['het']
            if results['het'][p]['exac_af']:
                # this patient has more than one var on this gene. copy it to hom_comp
                modes.append('hom_comp')
            populate_mode_p(results, 1, modes, p, exac_af, uclex_af, patient_db)
    return results

def populate_mode_p(results, copy,  modes, p, exac_af, uclex_af, patient_db):
    # copy: hom = both_het = 2, het = 1. 
    for mode in modes:
        results[mode][p] = results[mode].get(p, {'exac_af':[], 'uclex_af':[]})
        results[mode][p]['exac_af'].extend([exac_af]*copy)
        results[mode][p]['uclex_af'].extend([uclex_af]*copy)
        if 'hpo' not in results[mode]:
            # note that some patients don't have hpo features, and it will have 'All' as label. and some patients are not in the database, having label as None. useless to our purpose. remove this patient's entry!
            hpo = get_patient_observed_hpo(p, patient_db)
            # some of the patients hpo are tuples, not dics. report and convert!
            if not hpo or (len(hpo) == 1 and (not hpo[0][1] or hpo[0][1] == 'All')):
                del results[mode][p]
            else:
                results[mode][p]['hpo'] = hpo
            
'''
routes
'''
@app.route('/gene2/<gene_id>',methods=['GET'])
def gene_page2(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    hpo=request.args.get('hpo')
    hpo_freq = get_hpo_freq('hpo_freq.tsv')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    # get lof force
    # lof_p_hpo = get_lof_p_hpo(gene_id, db, patient_db)
    # rare_p_hpo = get_rare_var_p_hpo(gene_id, db, patient_db)
    # force_test = json.dumps(draw_force_graph(rare_p_hpo['het'], hpo_freq, hpo_db))
    inf = open('USH2A.force', 'r')
    force_test = json.load(inf)
    print force_test
    return render_template('gene2.html', 
            force_test = force_test)

# get sequence given region, and highlight the region. useful for design primers
@app.route('/sequence')
def sequence():
    var_id = request.args.get('variant_id')
    symbol = request.args.get('symbol')
    build = request.args.get('build') or 'grch37'
    (chrom,start,ref) = var_id.split('-',3)[:3]
    start = int(start)
    end = max(start + len(ref) - 1, start)
    strand = request.args.get('strand') or '1'
    expand_3prime = int(request.args.get('expand_3prime') or 450)
    expand_5prime = int(request.args.get('expand_5prime') or 450)

    # get sequence
    server = "http://%s.rest.ensembl.org" % build
    ext = '''/sequence/region/human/%(chrom)s:%(start)s..%(end)s:%(strand)s?expand_5prime=%(expand_5prime)s;expand_3prime=%(expand_3prime)s;''' % locals()
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
    region = [expand_5prime - 100, end - start + 200 ]
    seq = str(decoded['seq'])
    primer =  primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': 'phenopolis',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': region
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 59.0,
            'PRIMER_MIN_TM': 55.0,
            'PRIMER_MAX_TM': 67.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE':[[100,500],[400,1000]]
        })
    if 'PRIMER_RIGHT_0' not in primer:
        # this is not a very good return. will fix later
        return 'Cannot pick any primers for the given sequence'
    # formulate sequence
    seq = seq[:primer['PRIMER_RIGHT_0'][0]-primer['PRIMER_RIGHT_0'][1]+1] + '<span class="primer">' + seq[primer['PRIMER_RIGHT_0'][0]-primer['PRIMER_RIGHT_0'][1]+1:primer['PRIMER_RIGHT_0'][0]+1] + '</span>' + seq[primer['PRIMER_RIGHT_0'][0]+1:]
    seq = seq[:expand_5prime] + '<span class="highlight">' + seq[expand_5prime:expand_5prime+end-start+1] + '</span>' + seq[expand_5prime+end-start+1:]
    seq = seq[:primer['PRIMER_LEFT_0'][0]] + '<span class="primer">' + seq[primer['PRIMER_LEFT_0'][0]:primer['PRIMER_LEFT_0'][0]+primer['PRIMER_LEFT_0'][1]] + '</span>' + seq[primer['PRIMER_LEFT_0'][0]+primer['PRIMER_LEFT_0'][1]:]
    
    # construct object
    result = {'seq': seq,
            'var_id': var_id,
            'symbol': symbol,
            'left': primer['PRIMER_LEFT_0_SEQUENCE'],
            'left_tm': primer['PRIMER_LEFT_0_TM'],
            'left_gc': primer['PRIMER_LEFT_0_GC_PERCENT'], 
            'right': primer['PRIMER_RIGHT_0_SEQUENCE'],
            'right_tm': primer['PRIMER_RIGHT_0_TM'],
            'right_gc': primer['PRIMER_RIGHT_0_GC_PERCENT'],
            'product_length': primer['PRIMER_PAIR_0_PRODUCT_SIZE'],
            'expand_5prime': expand_5prime,
            'expand_3prime': expand_3prime
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

