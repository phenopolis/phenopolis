from __future__ import print_function
# Uncomment to run this module directly. TODO comment out.
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import runserver
from flask import Flask, current_app
from views import neo4j_driver
from views import my_patients
from views import session
import helper
import json



class MyPatientsPageTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        runserver.app.config['DB_NAME_USERS'] = 'test_users'
        self.app = runserver.app.test_client()
        #helper.create_neo4j_demo_user()
        helper.login(self.app)
        helper.my_patients_neo4j_data()


    def tearDown(self):
        self.app.get('/logout', follow_redirects=True)


    def test_my_patients_page(self): 
        page = self.app.get('/my_patients', follow_redirects=True)
        assert page.status_code == 200 # NB this test doesn't wait for the data to load.
        

    def commetn_in_test_get_individuals(self): #TODO LMTW
        app = Flask(__name__)
        with app.test_request_context():
            data = my_patients.get_individuals('demo')
            assert data['data'][1][0] == 'person1' 
            assert data['data'][1][1] == 'M' 
            assert data['data'][1][3] == 0.69 
            assert data['data'][1][6][0] == 'TTLL5' 
            assert data['data'][1][2][0]['data']['name'] == 'Visual impairment'
  
    def test_temp_get_individuals(self): 
        app = Flask(__name__)
        with app.test_request_context():
            data = my_patients.get_individuals('demo')
            assert data.status == '200 OK'
            parsed_json = json.loads(data.data)
            assert parsed_json['result'][0]['first_name'] == 'Jodi'

     

if __name__ == '__main__':
    unittest.main()

