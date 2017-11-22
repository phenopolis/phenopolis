# Uncomment to run this module directly. TODO comment out.
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import runserver
from flask import Flask, current_app
from views import neo4j_driver
from views import individual
import helper




class IndividualPageTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        runserver.app.config['DB_NAME_USERS'] = 'test_users'
        self.app = runserver.app.test_client()
        helper.create_neo4j_demo_user()
        helper.login(self.app)
        helper.my_patients_neo4j_data()


    def tearDown(self):
        self.app.get('/logout', follow_redirects=True)


    def test_venn_json_page(self): 
        page = self.app.get('/venn_json/person2', follow_redirects=True)
        assert page.status_code == 200 


    def test_get_feature_venn(self): 
        app = Flask(__name__)
        with app.test_request_context():
            feature_venn = individual.get_feature_venn('person2')
    

if __name__ == '__main__':
    unittest.main()

