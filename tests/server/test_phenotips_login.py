
import unittest

# Uncomment to run this module directly.
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
# End of Uncomment to run this module directly.

import runserver
import phenotips_python_client
from phenotips_python_client import PhenotipsClientNew
from config import config
import helper

from flask import Flask, session
from flask.ext.session import Session

import json


class PhenotipsLoginTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()
        helper.login(self.app)

    def tearDown(self):
        pass

    def test_login(self):
       
        if not config.LOCAL_WITH_PHENOTIPS:
            return
            
        conn = PhenotipsClientNew(test=True)
        phenotips_session = conn.request_phenotips_session('demo', 'demo123')
        assert(phenotips_session)

        invalid_user_session = conn.request_phenotips_session('demox', 'demo123')
        assert(not invalid_user_session)

        invalid_password_session = conn.request_phenotips_session('demo', 'demo123x')
        assert(not invalid_password_session)

        #auth='demo:demo123'
        #phenotips_session_auth = conn.request_phenotips_session(auth=auth)
        #assert(phenotips_session_auth)

        #username = (auth.split(':'))[0]
        #assert(username=='demo')
   
    def test_get_patient(self):
       
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)

        with self.app.session_transaction() as sess:

            conn.clear_cache()
            all_patients_from_phenotips = conn.get_patient(sess)
            assert(all_patients_from_phenotips)
            all_patients_from_phenotips = all_patients_from_phenotips.get('patientSummaries',[]) 
            all_eids = [p['eid'] for p in all_patients_from_phenotips if p['eid']]
            total_from_phenotips = len(all_eids)
            assert(total_from_phenotips>0)

            all_patients_from_cache = conn.get_patient(sess)
            assert(all_patients_from_cache)
            all_patients_from_cache = all_patients_from_cache.get('patientSummaries',[]) 
            all_eids = [p['eid'] for p in all_patients_from_cache if p['eid']]
            total_from_cache = len(all_eids)
            assert(total_from_cache == total_from_phenotips)

            new_phenotips_session = conn.request_phenotips_session('demo', 'demo123')
            sess['phenotips_session'] = new_phenotips_session 
            all_patients_new_session = conn.get_patient(sess)
            assert(all_patients_new_session)
            all_patients_new_session = all_patients_new_session.get('patientSummaries',[]) 
            all_eids = [p['eid'] for p in all_patients_new_session if p['eid']]
            total_new_session = len(all_eids)
            assert(total_from_cache == total_new_session)

            eid = 'P0000797'
            patient_from_phenotips = conn.get_patient(sess, eid)
            assert(patient_from_phenotips)
            eid_from_phenotips = str(patient_from_phenotips['external_id'])
            assert(eid_from_phenotips == eid)

            patient_from_cache = conn.get_patient(sess, eid)
            assert(patient_from_cache)
            eid_from_cache = str(patient_from_cache['external_id'])
            assert(eid_from_cache == eid_from_phenotips)

            authorised_patient_id = 'P0000002'
            permission = conn.get_permissions(sess, authorised_patient_id)
            assert(permission)

            unauthorised_patient_id = 'P0000001'
            permission = conn.get_permissions(sess, unauthorised_patient_id)
            assert(not permission)

    def test_get_patient_auth(self):
       
        ''' 
        As used by command line functions, e.g. load_patient
        '''
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)
        conn.clear_cache()
        auth='demo:demo123'
        session_phenotips = conn.create_session_with_phenotips(auth=auth)
        assert(session_phenotips)

        all_patients_from_phenotips = conn.get_patient(session_phenotips)
        assert(all_patients_from_phenotips)
        all_patients_from_phenotips = all_patients_from_phenotips.get('patientSummaries',[]) 
        all_eids = [p['eid'] for p in all_patients_from_phenotips if p['eid']]
        total_from_phenotips = len(all_eids)
        assert(total_from_phenotips>0)

        eid = 'P0000797'
        patient_from_phenotips = conn.get_patient(session_phenotips, eid)
        assert(patient_from_phenotips)
        eid_from_phenotips = str(patient_from_phenotips['external_id'])
        assert(eid_from_phenotips == eid)

    @staticmethod
    def load_patient(file_location):
        with open(file_location, 'r') as json_data:
            for line in json_data:
                dataset = json.loads(line)
            return dataset
     
    def test_update_patient(self):      
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)
        with self.app.session_transaction() as sess:
            file_location = "./tests/data/patient-update-name.json"
            patient = PhenotipsLoginTestCase.load_patient(file_location)
            eid = 'P0000005'
            conn.update_patient(eid, sess, patient)
            patient_from_phenotips = conn.get_patient(sess, eid)
            assert(patient_from_phenotips)
            patient_name = patient_from_phenotips['patient_name']
            name_from_phenotips = str(patient_name['first_name'])
            assert(name_from_phenotips == 'Updated')

            unauthorised_eid = 'P0000798'
            conn.update_patient(unauthorised_eid, sess, patient)            

    def test_create_and_delete_patient(self):
        if not config.LOCAL_WITH_PHENOTIPS:
            return
        conn = PhenotipsClientNew(test=True)
        with self.app.session_transaction() as sess:
            file_location = "./tests/data/simple-patient-Test01.json"
            patient = PhenotipsLoginTestCase.load_patient(file_location)
            conn.create_patient(sess, patient)
            eid = 'Test01'
            patient_retrieved = conn.get_patient(sess, eid)
            assert(patient_retrieved)
            patient_name = patient_retrieved['patient_name']
            name_retrieved = str(patient_name['first_name'])
            assert(name_retrieved == 'Paolo') 
            conn.delete_patient(eid, sess) 
            permissions = conn.get_permissions(session=session, eid=eid)
            assert(not permissions)
            
    def test_get_permissions(self):      
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)

        with self.app.session_transaction() as sess:  
            file_location = "./tests/data/simple-patient-Test01.json"
            patient = PhenotipsLoginTestCase.load_patient(file_location)
            conn.create_patient(sess, patient)
            eid = 'Test01'
            permissions = conn.get_permissions(sess, eid=eid)
            links = permissions['links']
            links_0 = links[0]
            allowed_methods = links_0['allowedMethods']
            assert(len(allowed_methods)==3)
            assert(str(allowed_methods[0]) == 'GET')
            assert(str(allowed_methods[1]) == 'PATCH')
            assert(str(allowed_methods[2]) == 'PUT')
            conn.delete_patient(eid, sess) 

            unauthorised_ID = 'P0000001'
            permissions = conn.get_permissions(sess, ID=unauthorised_ID)
            assert(not permissions)

    def test_update_permissions(self):      
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)

        with self.app.session_transaction() as sess:  
            file_location = "./tests/data/simple-patient-Test01.json"
            patient = PhenotipsLoginTestCase.load_patient(file_location)
            conn.create_patient(sess, patient)
            eid = 'Test01'
            original_permissions =  {"links":[{"allowedMethods":["GET","PATCH","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions","rel":"self"},{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/owner","rel":"https://phenotips.org/rel/owner"},{"allowedMethods":["DELETE","GET","PUT","PATCH"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/collaborators","rel":"https://phenotips.org/rel/collaborators"},{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/visibility","rel":"https://phenotips.org/rel/visibility"},{"allowedMethods":["DELETE","GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006","rel":"https://phenotips.org/rel/patientRecord"}],"owner":{"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/owner","rel":"https://phenotips.org/rel/owner"}],"id":"xwiki:XWiki.demo","name":"Demo Guest","email":"support@phenotips.org","type":"user"},"visibility":{"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/visibility","rel":"https://phenotips.org/rel/visibility"}],"level":"private","label":"private","description":"Private cases are only accessible to their owners, but they do contribute to aggregated statistics."},"collaborators":{"links":[{"allowedMethods":["DELETE","GET","PUT","PATCH"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/collaborators","rel":"https://phenotips.org/rel/collaborators"}],"collaborators":[]}}
            new_permissions =       {"links":[{"allowedMethods":["GET","PATCH","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions","rel":"self"},{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/owner","rel":"https://phenotips.org/rel/owner"},{"allowedMethods":["DELETE","GET","PUT","PATCH"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/collaborators","rel":"https://phenotips.org/rel/collaborators"},{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/visibility","rel":"https://phenotips.org/rel/visibility"},{"allowedMethods":["DELETE","GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006","rel":"https://phenotips.org/rel/patientRecord"}],"owner":{"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/owner","rel":"https://phenotips.org/rel/owner"}],"id":"xwiki:XWiki.demo","name":"Demo Guest","email":"support@phenotips.org","type":"user"},"visibility":{"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/visibility","rel":"https://phenotips.org/rel/visibility"}],"level":"public","label":"private","description":"Private cases are only accessible to their owners, but they do contribute to aggregated statistics."},"collaborators":{"links":[{"allowedMethods":["DELETE","GET","PUT","PATCH"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/collaborators","rel":"https://phenotips.org/rel/collaborators"}],"collaborators":[]}}
            conn.update_permissions(new_permissions, sess, eid=eid)
            permissions = conn.get_permissions(sess, eid=eid)
            assert(permissions['visibility']['label'] == 'public')
            conn.update_permissions(original_permissions, sess, eid=eid)
            permissions = conn.get_permissions(sess, eid=eid)
            assert(permissions['visibility']['label'] == 'private')
            conn.delete_patient(eid, sess) 

    def test_update_owner(self):      
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)

        with self.app.session_transaction() as sess: 
            file_location = "./tests/data/simple-patient-Test01.json"
            patient = PhenotipsLoginTestCase.load_patient(file_location)
            conn.create_patient(sess, patient)
            eid = 'Test01' 
            original_owner =    {"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions/owner","rel":"self"},{"allowedMethods":["GET","PATCH","PUT"],"href":"http://localhost:8080/rest/patients/P0000006/permissions","rel":"https://phenotips.org/rel/permissions"},{"allowedMethods":["DELETE","GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000006","rel":"https://phenotips.org/rel/patientRecord"}],"id":"xwiki:XWiki.demo","name":"Demo Guest","email":"support@phenotips.org","type":"user"}
            new_owner =         {"links":[{"allowedMethods":["GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000001/permissions/owner","rel":"self"},{"allowedMethods":["GET","PATCH","PUT"],"href":"http://localhost:8080/rest/patients/P0000001/permissions","rel":"https://phenotips.org/rel/permissions"},{"allowedMethods":["DELETE","GET","PUT"],"href":"http://localhost:8080/rest/patients/P0000001","rel":"https://phenotips.org/rel/patientRecord"}],"id":"xwiki:XWiki.Admin","name":"Administrator","email":"support@phenotips.org","type":"user"}
            conn.update_owner(new_owner, sess, eid=eid)
            permissions = conn.get_permissions(sess, eid=eid)
            assert(permissions['owner']['name'] == 'Administrator')
            conn.update_owner(original_owner, sess, eid=eid)
            permissions = conn.get_permissions(sess, eid=eid)
            assert(permissions['owner']['name'] == 'Demo Guest')
            conn.delete_patient(eid, sess) 

    def test_get_vocabularies(self):      
        if not config.LOCAL_WITH_PHENOTIPS:
            return   
        conn = PhenotipsClientNew(test=True)

        with self.app.session_transaction() as sess: 
            vocab = 'terms/HP:0000556'
            r = conn.get_vocabularies(session=sess,vocabulary=vocab)
            assert(str(r['name']) == 'Retinal dystrophy')

if __name__ == '__main__':
    unittest.main()
