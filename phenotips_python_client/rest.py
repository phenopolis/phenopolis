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

    def __init__(self, host='localhost', port='8080',debug=True,print_requests=True,test=False):
        self.site='%s:%s'%(host,port,)
        conn = pymongo.MongoClient(host='localhost', port=27017)
        if not test:
            self.db=conn.cache
        else:
            self.db=conn['test_cache']

    def request_phenotips_session(self, username=None, password=None):
        auth='%s:%s' % (username, password,)
        encoded_auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%encoded_auth, 'Accept':'application/json'}
        url='http://%s/rest/patients?start=%d&number=%d' % (self.site,0,1)
        s = requests.Session()
        response = s.get(url, headers=headers)
        if response :
            return s
        else:
            return None

    def create_session_with_phenotips(self, auth=None):
        encoded_auth=b2a_base64(auth).strip()
        headers={'Authorization':'Basic %s'%encoded_auth, 'Accept':'application/json'}
        url='http://%s/rest/patients?start=%d&number=%d' % (self.site,0,1)
        s = requests.Session()
        response = s.get(url, headers=headers)
        if response :
            username = (auth.split(':'))[0]
            session_dict = {'phenotips_session': s, 'user': username}
            return session_dict
        else:
            return None

    def get_phenotips_session(self, session):
        if not session or not 'phenotips_session' in session:
            return None
        phenotips_session = session['phenotips_session']
        if not phenotips_session:
            return None
        return phenotips_session

    def clear_cache(self):
        self.db.phenotips_cache.remove() 

    def get_patient(self,session,eid=None,number=10000,start=0):
        """
        Get patient with eid or all patients if not
        specified
        """
        phenotips_session = self.get_phenotips_session(session)
        if not phenotips_session:
            return None
        username = str((session['user']))

        s = phenotips_session

        headers={'Accept':'application/json'} 
        if not eid:
            url='http://%s/rest/patients?start=%d&number=%d' % (self.site,start,number)
            k={'url':url}
            k.update({'user':'%s'%username})
            k.update(headers)
            k = hashlib.md5(bencode.bencode(k)).hexdigest()
            r=self.db.phenotips_cache.find_one({'key':k})

            if r:
                return r
            else:
                r=s.get(url, headers=headers)
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
            k.update({'user':'%s'%username})
            k.update(headers)
            k = hashlib.md5(bencode.bencode(k)).hexdigest()
            r=self.db.phenotips_cache.find_one({'key':k})
            if r:
                return r
            else:
                r=s.get(url, headers=headers)
                try:
                    r=r.json()
                    r.update({'key':k})
                    self.db.phenotips_cache.insert_one(r)
                    return r
                except:
                    return None

    def patient_exists(self,session,eid):
        p=self.get_patient(session,eid)
        if p is None:
            return False
        else:
            return True

    def get_permissions(self,session,ID=None, eid=None):
        """
        Retrieves all permissions: owner, collaborators, visibility.
        """

        s = self.get_phenotips_session(session)
        if not s:
            return None
        if not ID:
            p=self.get_patient(session=session,eid=eid)
            ID=p['id']
        headers={'Accept':'application/json; application/xml'}
        r=s.get('http://%s/rest/patients/%s/permissions' % (self.site,ID), headers=headers)
        if not r:
            return None
        return r.json()

    # create patient
    def create_patient(self, session, patient):
        s = self.get_phenotips_session(session)
        if not s:
            return None
        headers={'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(patient,io)
        json_patient=io.getvalue()
        s.post('http://%s/rest/patients' % (self.site), headers=headers, data=json_patient)

    def update_patient(self, eid, session, patient):
        """
        Update patient if exists, otherwise create.
        """
        s = self.get_phenotips_session(session)
        if not s:
            return None
        patient['external_id']=eid
        if self.patient_exists(session=session,eid=eid):
            io=StringIO()
            json.dump(patient,io)
            json_patient=io.getvalue()
            print('update')
            print(json_patient)
            headers={'Content-Type':'application/json', 'Accept':'application/json'}
            s.put('http://%s/rest/patients/eid/%s' % (self.site,eid), headers=headers, data=json_patient)
        else:
            print('create')
            print(patient)
            self.create_patient(session=session,patient=patient)


    def update_permissions(self, permissions, session, ID=None, eid=None):
        """
        Update permissions of patient.
        """
        #permission = { "owner" : { "id" : "xwiki:XWiki.RachelGillespie" }, "visibility" : { "level":  "private" }, "collaborators" : [{ "id" : "xwiki:XWiki.UKIRDC", "level" : "edit" }, { "id" : "xwiki:Groups.UKIRDC Administrators)", "level" : "edit" }] }
        s = self.get_phenotips_session(session)
        if not s:
            return None
        if not ID:
            p=self.get_patient(session=session,eid=eid)
            ID=p['id']
        headers={'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(permissions,io)
        json_permissions=io.getvalue()
        p=s.put('http://%s/rest/patients/%s/permissions'% (self.site,ID), headers=headers, data=json_permissions, )
        print(p)
        return(p)


    def update_owner(self, owner, session, ID=None, eid=None):
        """
        Update owner of patient.
        """
        #permission = { "owner" : { "id" : "xwiki:XWiki.RachelGillespie" }, "visibility" : { "level":  "private" }, "collaborators" : [{ "id" : "xwiki:XWiki.UKIRDC", "level" : "edit" }, { "id" : "xwiki:Groups.UKIRDC Administrators)", "level" : "edit" }] }
        s = self.get_phenotips_session(session)
        if not s:
            return None
        if not ID:
            p=self.get_patient(session=session,eid=eid)
            ID=p['id']
        headers={'Content-Type':'application/json', 'Accept':'application/json'}
        io=StringIO()
        json.dump(owner,io)
        json_owner=io.getvalue()
        p=s.put('http://%s/rest/patients/%s/permissions/owner'%(self.site,ID), headers=headers, data=json_owner)
        print(p)
        return(p)


    def delete_patient(self, eid, session):
        """
        Delete patient.
        """
        s = self.get_phenotips_session(session)
        if not s:
            return None
        headers={'Content-Type':'application/json', 'Accept':'application/json'}
        p=s.delete('http://%s/rest/patients/eid/%s'%(self.site,eid), headers=headers)
        print(p)

    def update_phenotips_from_csv(self, info, auth, owner_group=[], collaborators=[], contact={}):
        """
        Each column in the csv file
        represent a patient atribute.
        This only supports one feature for now
        """
        info=pandas.read_csv(info,sep=',')
        print(info.columns.values)
        session = self.create_session_with_phenotips(auth=auth)
        for i, r, in info.iterrows():
            print(r)
            #if r['owner']!=owner: continue
            patient=dict()
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
            #update_patient(ID=r['sample'],session=session,patient=patient)
            self.update_patient(patient['external_id'], session, patient)
            #delete_patient(ID=r['sample'],session=session,patient=patient)
            # if patient exists, update patient, otherwise create patient
            #self.update_patient(eid=patient['external_id'],session=session,patient=patient)
            permissions = { "owner" : owner_group, "visibility" : { "level":  "private" }, "collaborators" : collaborators  }
            print(permissions)
            #self.update_permissions(permissions=permissions,eid=patient['external_id'],session=session)
            self.update_owner(owner=owner_group,session=session,eid=patient['external_id'])

    def patient_hpo(self, eid, session):
        """
        Retrieve HPO terms for patient
        """
        patient=self.get_patient(session=session,eid=eid)
        if patient:
            if 'features' in patient: return [f['id'] for f in patient['features']]
            else:  return []
        else: return []

    def dump_hpo_to_tsv(self, outFile, auth):
        """
        Dumps the HPO terms from a patient record
        to tsv file.
        """
        session = self.create_session_with_phenotips(auth=auth)
        patients=self.get_patient(session=session)['patientSummaries']
        #file(sprintf('uclex_hpo_%d-%d-%d.txt'),)
        hpo_file=open(outFile, 'w+')
        print('eid', 'hpo', 'genes', 'solved', sep='\t',file=hpo_file)
        for p in patients:
            eid=p['eid']
            print(eid)
            patient=self.get_patient(session=session,eid=eid)
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
        session = self.create_session_with_phenotips(auth=auth)
        patients=self.get_patient(session=session)['patientSummaries']
        for p in patients:
            eid=p['eid']
            print(eid)
            patient=self.get_patient(session=session,eid=eid)
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
        session = self.create_session_with_phenotips(auth=auth)
        patients=self.get_patient(session)['patientSummaries']
        for p in patients:
            eid=p['eid']
            print(eid)
            p=self.get_patient(session=session,eid=eid)
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
        session = self.create_session_with_phenotips(auth=auth)
        for eid in patient_ids:
            print(eid)
            p=self.get_patient(session=session,eid=eid)
            print(p)
            if p is None: raise 'patient does not exist maybe your credential are wrong?'
            # if patient does not exist in mongodb, create it
            if db.patients.find_one({'external_id':eid}) is None:
                db.patients.insert(p,w=0)
                continue
            for u in update_fields:
                print( u )
                if u not in p: p[u]=[]
                db.patients.update({'external_id':eid},{'$set':{u:p[u]}},w=0)


    def get_vocabularies(self,session,vocabulary):
        s = self.get_phenotips_session(session)
        if not s:
            return None
        # get vocabularies
        #http://localhost:1235/rest/vocabularies/terms/HP:0000556
        headers={'Accept':'application/json; application/xml'}
        r=s.get('http://%s/rest/vocabularies/%s'%(self.site,vocabulary), headers=headers)
        if not r:
            return None
        return r.json()




