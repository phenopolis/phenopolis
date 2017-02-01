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
        host = OFFLINE_CONFIG['mongodb']['DB_HOST'],
        port = int(OFFLINE_CONFIG['mongodb']['DB_PORT']),
    )
    if not test:
        return {
            'hpo_db': conn[OFFLINE_CONFIG['mongodb']['DB_NAME_HPO'],
            'phenopolis_db': conn[OFFLINE_CONFIG['mongodb']['DB_NAME']],
            'patient_db': conn[OFFLINE_CONFIG['mongodb']['DB_NAME_PATIENTS'],
        }

    else:
        return {
            'hpo_db': conn[test['hpo_db']],
            'phenopolis_db': conn[test['phenopolis_db']],
            'patient_db': conn[test['patient_db']],
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
