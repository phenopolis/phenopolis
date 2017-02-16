
import unittest


# TODO LMTW remove - for debugging
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import phenotips_python_client
# end TODO LMTW remove

import runserver
from phenotips_python_client import PhenotipsClient
from phenotips_python_client import PhenotipsClientNew
import httpie
from binascii import b2a_base64, a2b_base64
import requests

from flask import Flask, session
from flask.ext.session import Session
from flask import current_app


class PhenotipsLoginTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()

    def tearDown(self):
        pass

    def test_login(self):
       

        print('doing session')
        url = 'http://localhost:8080/rest/patients/P0000001/permissions'
        
        s = requests.Session()
        s.auth = requests.auth.HTTPDigestAuth('Admin', 'admin')
        r2 = s.get('http://localhost:8080/rest/patients/P0000001/permissions')


        encoded_auth=b2a_base64('Admin:admin').strip()
        headers={'Authorization':'Basic %s'%encoded_auth, 'Accept':'application/json'}
        r=requests.get(url, headers=headers)
        assert(r.status_code == 200) 

        s = requests.Session()
        r = s.get(url, headers=headers)
        assert(r.status_code == 200) 

        small_headers={'Accept':'application/json'}
        r = s.get(url, headers=small_headers)
        status_code = r.status_code 

        new_url='http://localhost:8080/rest/patients?start=0&number=5'
        r = s.get(new_url, headers=small_headers)
        status_code = r.status_code 


        with self.app.session_transaction() as sess:
            sess['phenotips'] = s 
            s2 = sess['phenotips']
            r = s2.get(new_url, headers=small_headers)
            status_code = r.status_code 
            assert(r.status_code == 200) 

        conn = PhenotipsClientNew()
        phenotips_session = conn.get_session('Admin', 'admin')
        with self.app.session_transaction() as sess:
            sess['phenotips'] = phenotips_session 
            s2 = sess['phenotips']
            r = s2.get(new_url, headers=small_headers)
            status_code = r.status_code 
            assert(r.status_code == 200) 

            conn.clear_cache()
            all_patients=conn.get_patient(sess['phenotips']).get('patientSummaries',[]) # TODO LMTW throw error if can't log in. ??
            all_eids=[p['eid'] for p in all_patients if p['eid']]
            total=len(all_eids)
            print('TOTAL NUMBER OF PATIENTS',total)

        temp =''


    
if __name__ == '__main__':
    unittest.main()
