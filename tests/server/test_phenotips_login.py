
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
        phenotips_session = conn.get_session('demo', 'demo123')
        assert(phenotips_session)

        invalid_user_session = conn.get_session('demox', 'demo123')
        assert(not invalid_user_session)

        invalid_password_session = conn.get_session('demo', 'demo123x')
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

            new_phenotips_session = conn.get_session('demo', 'demo123')
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


    
if __name__ == '__main__':
    unittest.main()
