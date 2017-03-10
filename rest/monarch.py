
import requests
import json
import pymongo
import sys

client = pymongo.MongoClient(port=27017)
patients_db = client['patients']


def score_phenotype(external_id):
    p=patients_db.patients.find_one({'external_id':external_id})
    x=dict()
    x['id']=p['report_id']
    x['features']=p['features']
    for f in p['features']:
        f['isPresent']={'yes':'true','no':'false'}[f['observed']]
        del f['observed']
        del f['label']
        del f['type']
    #print x
    #x={"features" : [{"id":"HP:0000505",'isPresent' : True},{"id":"HP:0000479", 'isPresent' : False}, {"id":"HP:0001010",'isPresent' : True},{"id":"HP:0000044", 'isPresent' : False}]}
    x=json.dumps(x)
    url='https://monarchinitiative.org/score/?annotation_profile={}'.format(x)
    print url
    r=requests.get(url,headers={'Content-Type':'application/json'})
    return r.json()

def hpo_terms(p):
    x=dict()
    x['id']=p['report_id']
    x['features']=p['features']
    for f in p['features']:
        f['isPresent']={'yes':'true','no':'false'}[f['observed']]
        del f['observed']
        del f['label']
        del f['type']
    return x['features']


def compare(individual,individual2):
    individual=conn.get_patient(eid=individual,session=session)
    individual2=conn.get_patient(eid=individual2,session=session)
    hpo1='+'.join([h['id'] for h in hpo_terms(individual)])
    hpo2='+'.join([h['id'] for h in hpo_terms(individual2)])
    url='https://monarchinitiative.org/compare/{}/{}.json'.format(hpo1,hpo2)
    print(url)
    r=requests.get(url,headers={'Content-Type':'application/json'})
    return r.json()

def get_phenotype_score(hpo):
    client = pymongo.MongoClient(port=27017)
    db = client['patients']
    for ind in db.patients.find():
        hpo2='+'.join([h['id'] for h in hpo_terms(ind)])
        url='https://monarchinitiative.org/compare/{}/{}.json'.format(hpo,hpo2)
        r=requests.get(url,headers={'Content-Type':'application/json'})
        x=r.json()
        if 'b' not in x:
            print 'SCORE', hpo, ind['external_id'], 0
        else:
            print 'SCORE', hpo, ind['external_id'], x['b'][0]['score']['score']



