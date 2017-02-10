'''
some common functions used across different analyses
'''
from __future__ import print_function, division
import pymongo
import ConfigParser
import os
import errno

'''
constants
'''
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X','Y']

'''
parse config file, and make config global. If test, set DB_HOST as 'localhost'
'''
def _parse_config():
    # return {'section':{'key1':'value1'...},...}
    config = ConfigParser.ConfigParser()
    
    # get this path, and then read common.cfg
    path = os.path.dirname(os.path.abspath(__file__))
    config.read(os.path.join(path,'common.cfg'))
    result = {}
    for section in config.sections():
        options = config.options(section)
        result[section] = {}
        for option in options:
            result[section][option] = config.get(section, option)
    return result

OFFLINE_CONFIG = _parse_config()
'''
get useful mongo collections
'''
def get_mongo_collections(test=None):
    conn = pymongo.MongoClient(
        host = OFFLINE_CONFIG['mongodb']['db_host'],
        port = int(OFFLINE_CONFIG['mongodb']['db_port']),
    )
    if not test:
        return {
            'hpo_db': conn[OFFLINE_CONFIG['mongodb']['db_name_hpo']],
            'phenopolis_db': conn[OFFLINE_CONFIG['mongodb']['db_name']],
            'patient_db': conn[OFFLINE_CONFIG['mongodb']['db_name_patients']],
            'pubmedbatch': conn[OFFLINE_CONFIG['mongodb']['db_name_pubmedbatch']],
        }

    else:
        return {
            'hpo_db': conn[test['hpo_db']],
            'phenopolis_db': conn[test['phenopolis_db']],
            'patient_db': conn[test['patient_db']],
            'pubmedbatch': conn[test['pubmedbatch_db']],
        }

'''
given chromosomes and db, return genes
'''
def get_chrom_genes(chroms, db):
    # give chrom numbers, get all genes on them
    result = []
    for chrom in chroms:
        chrom = str(chrom)
        if chrom not in VALID_CHROMOSOMES:
            raise ValueError('Error: %s is not a valid chromosome!' % chrom)
        genes = [g['gene_id'] for g in db.genes.find({'chrom':chrom})]
        result.extend(genes)
    return result

'''
mkdir -p
http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
basically, it only makes one system call, therefore avoid racing problems.
not required for python3, can use `os.makedirs(name,mode,exist_ok=True)
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

'''
get candidate genes and patients' hpos, solve, candidate genes, sex
'''
def get_candidate_genes(db, genes=None, fields=None):
    # set up some defaults. hpos = observed features.
    # solve would be 0 for unsolved and 1 for solved
    # sex 0 unknown, 1 male, 2 female
    SEX_DICT = {
            'F': 2,
            'M': 1,
            'U': 0,
        }
    SOLVE_DICT = {
            'solved':1,
            'unsolved':0,
        }

    fields = fields or ['hpo','solve','candidate_genes','sex']
    all_valid_p = [p for p in db.patients.find({}) if p.get('genes',[])]
    result = {}
    for p in all_valid_p:
        for g in p['genes']:
            # deal with hpo and solve and sex
            temp  = {f:p.get(f,None) for f in fields}
            if 'hpo' in fields:
                temp['hpo'] = [f for f in p['features'] if f['observed'] == 'yes']
            if 'solve' in fields:
                temp['solve'] = SOLVE_DICT[p['solved']['status']]
            if 'sex' in fields:
                temp['sex'] = SEX_DICT[p['sex']]

            result[g['gene']] = result.get(g['gene'],[])
            result[g['gene']].append(temp)
    return result
