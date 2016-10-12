#!/bin/env python
'''
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

cutoff freq on common ancestor: 0.1
the node with smallest freq above 0.1 will be set as a group
When a patient belongs to different groups, colouring can be difficult
    we now try to assign it to the group where this patient's minimal hpo are the most
    if tie, random assign to tied group
optional to nearest and average
'''

from optparse import OptionParser
from lookups import *
import pymongo
import json
import math

conn = pymongo.MongoClient(host='localhost', port=27017)
hpo_db=conn['hpo']
db = conn['uclex-old']
patient_db=conn['patients']
###############################
# parse options
###############################
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-c", "--cutoff",
                  dest="cutoff", default='0.1',
                  help="ignore connections with values bigger than or equal to cutoff. this value is also used to colour the groups [default: %default]")
parser.add_option("-m", "--mode", dest="mode", default='nearest',
                  help='either nearest or average [default: %default]')
                  
(options, args) = parser.parse_args()
'''
defs
'''
def draw_force_graph(patients, hpo_freq, hpo_db):
    '''
    patients = {p_id:{'hpo':[(HP:001, hel), ...]}}
    result = '{
      "nodes":[
          {"name":"Myriel","group":1, "radius": 5, "hpo":[["HP:0000345","hell"],["HP..]] },
          {"name":"Napoleon","group":1, "radius": 5, "hpo": ..},
          {"name":"Mlle.Baptistine","group":1, "radius": 5, "hpo": ..},
          {"name":"Mme.Magloire","group":1, "radius": 5, "hpo": ..}
      ],
      "links":[{"source":1,"target":0,"value":1},
      {"source":2,"target":0,"value":8},
      {"source":3,"target":0,"value":10},
      ]
  }
    '''
    hpo_groups = {}
    # convert patients to array to conserve order, and populate hpo_groups
    ii = 1 # for debug. can delete safely.
    for p in patients:
        ii += 1
        # replace obsolete hpo first
        temp_hpos = [replace_hpo(hpo_db, i) for i in patients[p]['hpo']]
        # minimise
        min_set = hpo_minimum_set(hpo_db, hpo_ids=[i[0] for i in temp_hpos])
        patients[p]['hpo'] = [list(i) for i in temp_hpos if i[0] in min_set]
        hpo_groups = hpo_grouping(hpo_groups, {'p_id':p, 'hpo':patients[p]['hpo']}, hpo_freq, hpo_db, float(options.cutoff))

    patients = [{'name':key, 'group':0, 'hpo':patients[key]['hpo']} for key in patients]
    result = {'nodes':patients, 'links':[]}
    p_combo = itertools.combinations(range(len(patients)), 2)
    for pair in p_combo:
        p1 = [i[0] for i in patients[pair[0]]['hpo']]
        p2 = [i[0] for i in patients[pair[1]]['hpo']]
        if options.mode == 'nearest':
            value = patients_hpo_nearest_similarity(p1, p2, hpo_freq, hpo_db, cutoff=float(options.cutoff))
        else:
            value = patients_hpo_average_similarity(p1, p2, hpo_freq, hpo_db, cutoff=float(options.cutoff))
        if not value:
            continue
        result['links'].append({
            'source':pair[0],
            'target':pair[1],
            'value':value
            })
    check_hpo_grouping(hpo_groups, hpo_db, hpo_freq)
    print hpo_groups
    update_result(result, hpo_groups, hpo_freq, hpo_db)

    return result

def hpo_grouping(local_hpo_groups, patient, hpo_freq, hpo_db, cutoff):
    # group patient to hpo groups, and update statistics
    '''
    hpo_groups = {'HP:0000234': {'group':0, 'patients': [p_ids]}}
    patient = {'p_id':'foo', 'hpo':[(HP:001, hel), ...]} should already have been minimised
    '''
    for hpo in patient['hpo']:
        # get ancestors
        ancestors = get_hpo_ancestors_array(hpo_db, hpo[0])
        # some ancestors are already in the hpo_groups?
        chosen  = set(ancestors) & set(local_hpo_groups.keys())
        if chosen:
            for c in chosen:
                local_hpo_groups[c]['patients'].append(patient['p_id'])
        else:
            # not in the hpo_group, need to find the left-closest to cutoff
            the_min = 0
            the_one = ''
            the_max = 0
            the_min_one = ''
            for a in ancestors:
                if the_max < hpo_freq[a]['freq'] <= cutoff:
                    the_one = a
                    the_max = hpo_freq[a]['freq']
                if hpo_freq[a]['freq'] > the_min:
                    the_min = hpo_freq[a]['freq']
                    the_min_one = a
            if not the_one:
                # big group, none of the ancestors below the cutoff
                # find the ancestor with the min freq
                the_one = the_min_one

            local_hpo_groups[the_one] = {'patients': [patient['p_id']]}
    return local_hpo_groups

def check_hpo_grouping(hpo_groups, hpo_db, hpo_freq):
    # this make sure group HPOs are independent of each other
    # some groups might be ancestors of others. merge them to the ancestors,
    #  and calculate grouping
    bad = []
    total = []
    for i in hpo_groups:
        total.extend(hpo_groups[i]['patients'])
        if not i:
            print 'hell'
            continue
        ancestors = get_hpo_ancestors_array(hpo_db, i)
        for j,v in hpo_groups.iteritems():
            if i == j:
                continue
            if j in ancestors:
                bad.append(i)
                v['patients'].extend(hpo_groups[i]['patients'])
    for b in bad:
        if b not in hpo_groups:
            print 'odd, %s already deleted from hpo_groups?' % b
            continue
        del hpo_groups[b]
    # set group number, remove repetitive patients, and calculate total appearance
    # set radius: r = sqrt(observed_patients_frequency / expected_patients_frequency) * 5 
    total = float(len(set(total)))
    group = 0
    for i,v in hpo_groups.iteritems():
        group += 1
        v['patients'] = set(v['patients'])
        v['group'] = group
        observed = len(v['patients']) / total
        print i
        print observed
        print hpo_freq[i]['freq']
        v['radius'] = math.sqrt(observed / hpo_freq[i]['freq']) * 5
        print v['radius']


def update_result(result, hpo_groups, hpo_freq, hpo_db):
    # update the result with hpo grouping, 
    '''
    result = '{
      "nodes":[
          {"name":"Myriel","group":1, "radius": 5, "hpo":[["HP:0000345","hell"],["HP..]] },
          {"name":"Napoleon","group":1, "radius": 5, "hpo": ..},
          {"name":"Mlle.Baptistine","group":1, "radius": 5, "hpo": ..},
          {"name":"Mme.Magloire","group":1, "radius": 5, "hpo": ..}
      ],
      "links":[{"source":1,"target":0,"value":1},
      {"source":2,"target":0,"value":8},
      {"source":3,"target":0,"value":10},
      ]
  }
    '''
    group_hpos = set(hpo_groups.keys())
    for v in result['nodes']:
        # see which group it belongs
        groups = {} # {'HP:0000345': 2, ...}
        maxi = 0
        the_groups = []
        for h in v['hpo']:
            # get the groups. note that one hpo can belong to two or more groups
            #  that are independent with each other
            this_groups = set(get_hpo_ancestors_array(hpo_db, h[0])) & group_hpos
            for g in this_groups:
                groups[g] = groups.get(g, 0) + 1
                if groups[g] > maxi:
                    the_groups = [g]
                elif groups[g] == maxi:
                    the_groups.append(g)
        # random select from the_groups to assign
        temp = random.choice(the_groups)
        v['group'] = hpo_groups[temp]['group']
        v['radius'] = hpo_groups[temp]['radius']
        
            

def patients_hpo_nearest_similarity(p1, p2, hpo_freq, hpo_db, cutoff=0.5):
    # calculate patients similarity given hpo terms
    # p1 = [HPid, HPid]
    '''
    hpo similarity is calculated by Lin's similarity
    find the most similar hpo terms between p1 and p2, and report the
        similarity score
    if nearest common ancestor's freq >= cutoff, no link
    '''
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
    t1 = hpo_freq[h1]['freq'] 
    t2 = hpo_freq[h2]['freq']
    # nearest ancestor, and get freq
    n_a = get_hpo_nearest_common_ancestors(hpo_db, h1, h2, hpo_freq)
    if float(hpo_freq[n_a[0]]['freq']) >= cutoff:
        return 0
    max_an_IC = IC(hpo_freq[n_a[0]]['freq'])
    return max_an_IC * 2 / (IC(t1) + IC(t2))

def get_rare_var_p_hpo(gene_id, db, patient_db):
    #return {'hom_comp':{p_id1:{hpo: [(HP:1234,hell)], exac_af:[0.0021], uclex_af:[0.001]},
    #        'het':{p_id2: {hpo:[(HP:2345,yeah)], exac_af:[0.001,0.002],uclex_af:[0.002,0.001]}}

    # sometimes variant is not in vcf. move it to debug/bad_variants for inspection and later clean
    # for het, cut at exac_af = 0.001
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
        this = vcf_query(variant_str=v)
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
    # clean het with cut at exac_af==0.001
    print 'before clean, rare_het length = %s' % len(results['het'])
    bad_p = []
    for p, v in results['het'].iteritems():
        if min(v['exac_af']) > 0.001:
            bad_p.append(p)
    for p in bad_p:
        del results['het'][p]
    print 'after clean, rare_het length = %s' % len(results['het'])
    print 'rare_hom_comp length = %s' % len(results['hom_comp'])
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
main
'''
# for the time being, experiment on USH2A. write to USH2A.force
if options.mode not in ['nearest', 'average']:
    raise('mode has to be either nearest or average')
genes = ['OPN4']
hpo_freq = get_hpo_size_freq('hpo_freq.tsv')
for gene_id in genes:
    outf = open(gene_id + '.force','w')
    if not gene_id.startswith('ENSG'): gene_id = get_gene_by_name(db, gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    rare_p_hpo = get_rare_var_p_hpo(gene_id, db, patient_db)
    force_test = json.dumps(draw_force_graph(rare_p_hpo['het'], hpo_freq, hpo_db))
    json.dump(force_test, outf)
