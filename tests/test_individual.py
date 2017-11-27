# Uncomment to run this module directly. TODO comment out.
#import sys, os
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import runserver
from flask import Flask, current_app
from views import neo4j_driver
from views import individual
import helper
import json




class IndividualPageTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        runserver.app.config['DB_NAME_USERS'] = 'test_users'
        self.app = runserver.app.test_client()
        helper.login(self.app)
        helper.my_patients_neo4j_data()


    def tearDown(self):
        self.app.get('/logout', follow_redirects=True)


    def test_venn_json_page(self): 
        page = self.app.get('/venn_json/person2', follow_redirects=True)
        assert page.status_code == 200 
        parsed_json = json.loads(page.data)
        # Check some key points
        i=4
        assert parsed_json['result'][i]['key'][0] == 'Abnormality of the retina'
        assert parsed_json['result'][i]['key'][1] == 'Visual impairment'
        assert parsed_json['result'][i]['value'][0] == 'TTLL5'
        assert parsed_json['result'][i]['value'][1] == 'DRAM2'
        i=5
        assert parsed_json['result'][i]['key'][0] == 'Macular dystrophy'
        assert parsed_json['result'][i]['key'][1] == 'Visual impairment'
        assert parsed_json['result'][i]['value'][0] == 'TRIM32'
        i=6
        assert parsed_json['result'][i]['key'][0] == 'Abnormality of the retina'
        assert parsed_json['result'][i]['key'][1] == 'Macular dystrophy'
        assert parsed_json['result'][i]['key'][2] == 'Visual impairment'
        assert not parsed_json['result'][i]['value']
   
if __name__ == '__main__':
    unittest.main()

