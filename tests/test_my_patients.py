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


    def test_my_patients_page(self): # TODO LMTW comment in
        page = self.app.get('/my_patients', follow_redirects=True)
        assert page.status_code == 200
        # Need to wait for data before it can be checked. TODO LMTW
        #assert 'person1' in page.data 
        #assert 'M' in page.data 
        #assert '0.69' in page.data 
        #assert 'Visual impairment' in page.data
        #assert 'TTLL5' in page.data

    def test_get_individuals(self):
        #app = Flask(__name__)
        #with app.test_request_context():
        #    with self.app.session_transaction() as session:
        #        session['user']='demo'
        #        #print('#####################')
        #        #print(temp)

        app = Flask(__name__)
        with app.test_request_context():
            data = my_patients.get_individuals()
            assert 'person1' in data 
            assert 'M' in data 
            assert '0.69' in data 
            assert 'Visual impairment' in data
            assert 'TTLL5' in data

if __name__ == '__main__':
    unittest.main()

