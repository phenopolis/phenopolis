from __future__ import print_function
import requests
from binascii import b2a_base64, a2b_base64
import re
import os.path
import sys 
import httplib
import urllib
import urlparse
import re
import socket
import gzip
import StringIO
from StringIO import StringIO
import time
import random
from binascii import b2a_base64, a2b_base64
#from email.utils import fix_eols
import json
import pandas
import hashlib
import bencode


#import requests_cache
#requests_cache.install_cache('phenotips_cache')
#from cachecontrol import CacheControl
#from cachecontrol.caches import FileCache
#sess = CacheControl(requests.Session(), cache=FileCache('.webcache'))
#sess = cachecontrol.CacheControl(requests.Session())

import pymongo
#requests_cache.backends.mongo.MongoCache(db_name='phenotips_cache', connectoin)

class PhenotipsClient():

    def __init__(self, host='localhost', port='8080',debug=True,print_requests=True):
        self.site='%s:%s'%(host,port,)
        conn = pymongo.MongoClient(host='localhost', port=27017)
        self.db=conn.cache

    def get_patient(self,auth,eid=None,number=10000,start=0,cache=True):
        """
        Get patient with eid or all patients if not
        specified
        """
        auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json'}
        if not eid:
            url='http://%s/rest/patients?start=%d&number=%d' % (self.site,start,number)
            k={'url':url}
            k.update(headers)
            k = hashlib.md5(bencode.bencode(k)).hexdigest()
            r=self.db.phenotips_cache.find_one({'key':k})
            if r and cache:
                print('CACHED')
                return r
            else:
                print('NOT CACHED')
                r=requests.get(url, headers=headers)
                try:
                    r=r.json()
                    r.update({'key':k})
                    self.db.phenotips_cache.insert_one(r)
                    return r
                except:
                    return None
        else:
            url='http://%s/rest/patients/eid/%s' % (self.site,str(eid))
            k={'url':url}
            k.update(headers)
            k = hashlib.md5(bencode.bencode(k)).hexdigest()
            r=self.db.phenotips_cache.find_one({'key':k})
            if r and cache:
                print('CACHED')
                return r
            else:
                print('NOT CACHED')
                r=requests.get(url, headers=headers)
                try:
                    r=r.json()
                    r.update({'key':k})
                    self.db.phenotips_cache.insert_one(r)
                    del r['_id']
                    return r
                except:
                    return None

    def patient_exists(self,auth,eid):
        p=self.get_patient(auth,eid)
        if p is None:
            return False
        else:
            return True

    def get_permissions(self,auth,ID):
        """
        Retrieves all permissions: owner, collaborators, visibility.
        """
        auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json; application/xml'}
        #p=self.get_page('/patients/%s/permissions',ID, headers=headers)
        r=requests.get('http://%s/rest/patients/%s/permissions' % (self.site,ID), headers=headers)
        return r.json()

    # create patient
    def create_patient(self, auth, patient):
        headers={'Authorization':'Basic %s'% b2a_base64(auth).strip(), 'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(patient,io)
        json_patient=io.getvalue()
        #p=self.get_page('/rest/patients', headers=headers, post=json_patient)
        print(p)
        return(p)

    def update_patient(self, eid, auth, patient):
        """
        Update patient if exists, otherwise create.
        """
        patient['external_id']=eid
        if self.patient_exists(auth=auth,eid=eid):
            io=StringIO()
            json.dump(patient,io)
            json_patient=io.getvalue()
            print('update')
            print(json_patient)
            headers={'Authorization':'Basic %s'% b2a_base64(auth).strip(),'Content-Type':'application/json', 'Accept':'application/json'}
            self.get_page('/rest/patients/eid/%s'%eid, headers=headers, post=json_patient, special='PUT')
        else:
            print('create')
            print(patient)
            self.create_patient(auth=auth,patient=patient)


    def update_permissions(self, permissions, auth, ID=None, eid=None):
        """
        Update permissions of patient.
        """
        #permission = { "owner" : { "id" : "xwiki:XWiki.RachelGillespie" }, "visibility" : { "level":  "private" }, "collaborators" : [{ "id" : "xwiki:XWiki.UKIRDC", "level" : "edit" }, { "id" : "xwiki:Groups.UKIRDC Administrators)", "level" : "edit" }] }
        if not ID:
            p=self.get_patient(auth=auth,eid=eid)
            ID=p['id']
        auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%auth, 'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(permissions,io)
        json_permissions=io.getvalue()
        p=self.get_page('/patients/%s/permissions'%ID, headers=headers, post=json_permissions, special='PUT')
        print(p)
        return(p)


    def update_owner(self, owner, auth, ID=None, eid=None):
        """
        Update owner of patient.
        """
        #permission = { "owner" : { "id" : "xwiki:XWiki.RachelGillespie" }, "visibility" : { "level":  "private" }, "collaborators" : [{ "id" : "xwiki:XWiki.UKIRDC", "level" : "edit" }, { "id" : "xwiki:Groups.UKIRDC Administrators)", "level" : "edit" }] }
        if not ID:
            p=self.get_patient(auth=auth,eid=eid)
            ID=p['id']
        auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%auth, 'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(owner,io)
        json_owner=io.getvalue()
        p=self.get_page('/patients/%s/permissions/owner'%ID, headers=headers, post=json_owner, special='PUT')
        print(p)
        return(p)


    def delete_patient(self, eid, auth):
        """
        Delete patient.
        """
        auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%auth, 'Content-Type':'application/json', 'Accept':'application/json'}
        p=self.get_page('/rest/patients/eid/%s'%eid, headers=headers, post='', special='DELETE')
        print(p)

    def update_phenotips_from_csv(self, info, auth, owner_group=[], collaborators=[], contact={}):
        """
        Each column in the csv file
        represent a patient atribute.
        This only supports one feature for now
        """
        info=pandas.read_csv(info,sep=',')
        print(info.columns.values)
        for i, r, in info.iterrows():
            print(r)
            #if r['owner']!=owner: continue
            patient=dict()
            #auth=login
            #auth=b2a_base64(auth).strip()
            patient['external_id']=r['sample']
            if  not isinstance(r['sample'],basestring) or len(r['sample']) < 4: continue
            if 'ethnicity' in r:
                ethnicity=r['ethnicity']
                if isinstance(ethnicity,str):
                    patient["ethnicity"]={"maternal_ethnicity":[ethnicity],"paternal_ethnicity":[ethnicity]}
                else:
                    patient["ethnicity"]={"maternal_ethnicity":[],"paternal_ethnicity":[]}
            patient["prenatal_perinatal_history"]={}
            patient["prenatal_perinatal_phenotype"]={"prenatal_phenotype":[],"negative_prenatal_phenotype":[]}
            patient['reporter']=r['owner']
            if 'gender' in r:
                gender=dict({'M':'M','m':'M','f':'F','F':'F','1':'M','2':'F'}).get(str(r['gender']),'U')
                patient['sex']=gender
            #if 'solved' in r: patient['solved']=r['solved']
            #patient['contact']={ "user_id":r['owner'], "name":r['owner'], "email":'', "institution":'' }
            patient['contact']=contact
            patient['clinicalStatus']={ "clinicalStatus":"affected" }
            #patient['disorders']=[ { "id":r['phenotype'], 'label':''} ]
            print(patient)
            r['phenotype']=str(r['phenotype'])
            patient['features']=[ { "id":hpo, 'label':'', 'type':'phenotype', 'observed':'yes' } for hpo in r['phenotype'].split(';') ]
            #update_patient(ID=r['sample'],auth=auth,patient=patient)
            self.update_patient(patient['external_id'], auth, patient)
            #delete_patient(ID=r['sample'],auth=auth,patient=patient)
            # if patient exists, update patient, otherwise create patient
            #self.update_patient(eid=patient['external_id'],auth=auth,patient=patient)
            permissions = { "owner" : owner_group, "visibility" : { "level":  "private" }, "collaborators" : collaborators  }
            print(permissions)
            #self.update_permissions(permissions=permissions,eid=patient['external_id'],auth=auth)
            self.update_owner(owner=owner_group,auth=auth,eid=patient['external_id'])

    def patient_hpo(self, eid, auth):
        """
        Retrieve HPO terms for patient
        """
        patient=self.get_patient(auth,eid=eid)
        if patient:
            if 'features' in patient: return [f['id'] for f in patient['features']]
            else:  return []
        else: return []

    def dump_hpo_to_tsv(self, outFile, auth):
        """
        Dumps the HPO terms from a patient record
        to tsv file.
        """
        patients=self.get_patient(auth)['patientSummaries']
        #file(sprintf('uclex_hpo_%d-%d-%d.txt'),)
        hpo_file=open(outFile, 'w+')
        print('eid', 'hpo', 'genes', 'solved', sep='\t',file=hpo_file)
        for p in patients:
            eid=p['eid']
            print(eid)
            patient=self.get_patient(auth,eid)
            print(patient)
            if 'features' in patient:
                hpo=','.join([f['id'] for f in patient['features']])
            else:
                hpo=''
            if 'genes' in patient:
                genes=','.join([g['gene'] for g in patient['genes']])
            else:
                genes=''
            if 'solved' in patient:
                solved=patient['solved']['status']
            else:
                solved='unsolved'
            print(eid, hpo, genes, solved, sep='\t',file=hpo_file)


    def dump_patient_to_json(self, auth):
        """
        Dumps patient to JSON.
        """
        auth='%s:%s' % (owner, password,)
        patients=self.get_patient(auth)['patientSummaries']
        for p in patients:
            eid=p['eid']
            print(eid)
            patient=self.get_patient(auth,eid)
            io=StringIO()
            json.dump(patient,io)
            json_patient=io.getvalue()
            print(json_patient)


    def dump_to_mongodb(self, auth, mongo_host='localhost', mongo_port=27016, mongo_dbname='patients'):
        """
        Dumps all patients to a Mongo database
        """
        import pymongo
        client = pymongo.MongoClient(host=mongo_host, port=int(mongo_port))
        db=client[mongo_dbname]
        db.patients.drop()
        patients=self.get_patient(auth)['patientSummaries']
        for p in patients:
            eid=p['eid']
            print(eid)
            p=self.get_patient(auth,eid)
            db.patients.insert(p,w=0)
        db.patients.ensure_index('external_id')
        db.patients.ensure_index('report_id')
        db.patients.ensure_index('features.id')
        db.patients.ensure_index('sex')
        db.patients.ensure_index('genes.gene')
        db.patients.ensure_index('solved')
        db.patients.ensure_index('clinicalStatus.clinicalStatus')
        db.patients.ensure_index('specificity.score')


    def update_mongodb(self, auth, mongo_host='localhost', mongo_port=27016, mongo_dbname='patients',update_fields=['features'],patient_ids=[]):
        """
        Will only update fields which have changed in Phenotips
        """
        import pymongo
        client = pymongo.MongoClient(host=mongo_host, port=int(mongo_port))
        db=client[mongo_dbname]
        for eid in patient_ids:
            print(eid)
            p=self.get_patient(auth,eid)
            print(p)
            if p is None: raise 'patient does not exists maybe your credential are wrongs?'
            # if patient does not exist in mongodb, create it
            if db.patients.find_one({'external_id':eid}) is None:
                db.patients.insert(p,w=0)
                continue
            for u in update_fields:
                print( u )
                if u not in p: p[u]=[]
                db.patients.update({'external_id':eid},{'$set':{u:p[u]}},w=0)


    def get_vocabularies(self,auth,vocabulary):
        auth=b2a_base64(auth).strip()
        # get vocabularies
        #http://localhost:1235/rest/vocabularies/terms/HP:0000556
        headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json; application/xml'}
        p=self.get_page('/rest/vocabularies/%s'%vocabulary, headers=headers)
        print(p)
        return p




