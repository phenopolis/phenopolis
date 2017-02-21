
import unittest

# Uncomment to run this module directly.
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
# End of Uncomment to run this module directly.

import runserver
import phenotips_python_client
from phenotips_python_client import PhenotipsClientNew
from phenotips_python_client import PhenotipsClient # TODO LMTW remove
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

            patient_id = 'P0000797'
            patient_from_phenotips = conn.get_patient(sess, patient_id)
            assert(patient_from_phenotips)
            eid_from_phenotips = str(patient_from_phenotips['external_id'])
            assert(eid_from_phenotips == patient_id)

            patient_from_cache = conn.get_patient(sess, patient_id)
            assert(patient_from_cache)
            eid_from_cache = str(patient_from_cache['external_id'])
            assert(eid_from_cache == eid_from_phenotips)

            authorised_patient_id = 'P0000002'
            permission = conn.get_permissions(sess, authorised_patient_id)
            assert(permission)

            unauthorised_patient_id = 'P0000001'
            permission = conn.get_permissions(sess, unauthorised_patient_id)
            assert(not permission)

    @staticmethod
    def single_patient():
        return {"genes": [{"status": "", "gene": "", "comments": ""}], "external_id": "P0000101", "features": [{"observed": "no", "type": "phenotype", "id": "HP:0000593", "label": "Abnormality of the anterior chamber"}, {"observed": "no", "type": "phenotype", "id": "HP:0000481", "label": "Abnormality of the cornea"}, {"observed": "yes", "id": "HP:0000479", "type": "phenotype", "qualifiers": [{"type": "age_of_onset", "id": "HP:0003577", "label": "Congenital onset"}, {"type": "laterality", "id": "HP:0012832", "label": "Bilateral"}], "label": "Abnormality of the retina"}, {"observed": "no", "type": "phenotype", "id": "HP:0000589", "label": "Coloboma"}, {"observed": "no", "type": "phenotype", "id": "HP:0000316", "label": "Hypertelorism"}, {"observed": "no", "type": "phenotype", "id": "HP:0000601", "label": "Hypotelorism"}, {"observed": "no", "type": "phenotype", "id": "HP:0000568", "label": "Microphthalmos"}, {"observed": "yes", "id": "HP:0000556", "type": "phenotype", "qualifiers": [{"type": "age_of_onset", "id": "HP:0003577", "label": "Congenital onset"}, {"type": "laterality", "id": "HP:0012832", "label": "Bilateral"}, {"type": "severity", "id": "HP:0012828", "label": "Severe"}], "label": "Retinal dystrophy"}, {"observed": "yes", "id": "HP:0000550", "type": "phenotype", "qualifiers": [{"type": "age_of_onset", "id": "HP:0003593", "label": "Infantile onset"}, {"type": "laterality", "id": "HP:0012832", "label": "Bilateral"}], "label": "Undetectable electroretinogram"}, {"observed": "yes", "id": "HP:0000505", "type": "phenotype", "qualifiers": [{"type": "age_of_onset", "id": "HP:0003577", "label": "Congenital onset"}, {"type": "laterality", "id": "HP:0012832", "label": "Bilateral"}, {"type": "severity", "id": "HP:0012828", "label": "Severe"}], "label": "Visual impairment"}]}
    
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
            eid = 'P0000003'
            conn.update_patient(eid, sess, patient)

    #def test_create_patient(self):
    #    if not config.LOCAL_WITH_PHENOTIPS:
    #        return
    #    conn = PhenotipsClientNew(test=True)
    #    with self.app.session_transaction() as sess:
    #        file_location = "./tests/data/simple-patient-P0000006.json"
    #        patient = PhenotipsLoginTestCase.load_patient(file_location)
    #        conn.create_patient(sess, patient)                 

if __name__ == '__main__':
    unittest.main()
