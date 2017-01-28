'''
some common functions used across different analyses
'''
from __future__ import print_function, division
import pymongo
import ConfigParser
import os

'''
parse config file, and make config global. If test, set DB_HOST as 'localhost'
'''
def _parse_config(test=None):
    # basic settings are under the `default` section
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
    if test:
        result['mongodb']['db_host'] = 'localhost'
    return result

OFFLINE_CONFIG = _parse_config()
'''
get useful mongo collections
'''
def get_mongo_collections():
    conn = pymongo.MongoClient(
            host = OFFLINE_CONFIG['mongodb']['db_host'],
            port = int(OFFLINE_CONFIG['mongodb']['db_port']),
            )
    return {
            'hpo_db': conn['hpo'],
            'phenopolis_db': conn[OFFLINE_CONFIG['mongodb']['db_name']],
            'patient_db': conn['patients'],
    }
