# Uncomment to run this module directly. TODO comment out.
#import sys, os
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import unittest
import subprocess
import runserver
from flask import Flask, current_app, jsonify
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

  
    def test_my_patients_functionality(self): 
        app = Flask(__name__)
        with app.test_request_context():
            records = my_patients.get_individuals('demo')
            # Here we create the Flask Response object, containing json, 
            # that the /my_patients page receives. We then test 
            # that the expected data is available.
            data=jsonify(result=records)
            assert data.status == '200 OK'
            parsed_json = json.loads(data.data)
            # First person.
            i=0
            assert parsed_json['result'][i]['individual'] == 'person2'
            assert parsed_json['result'][i]['gender'] == 'F'
            for pheno in parsed_json['result'][i]['phenotypes'] :
                assert (pheno['name'] == 'Abnormality of the retina' or 
                        pheno['name'] == 'Visual impairment' or 
                        pheno['name'] == 'Macular dystrophy')            
            assert parsed_json['result'][i]['phenotypeScore'] == 0.69
            assert parsed_json['result'][i]['hom_count'] == 1
            assert parsed_json['result'][i]['het_count'] == 2
            for gene in parsed_json['result'][i]['genes'] :
                assert gene == 'RPGR' or gene == 'TTLL5' or gene == 'DRAM2' or gene == 'TRIM32'
            # Next person.
            i=1
            assert parsed_json['result'][i]['individual'] == 'person1'
            assert parsed_json['result'][i]['gender'] == 'M'
            assert parsed_json['result'][i]['phenotypes'][0]['name'] == 'Visual impairment'
            assert parsed_json['result'][i]['phenotypeScore'] == 0.69
            assert parsed_json['result'][i]['hom_count'] == 1
            assert parsed_json['result'][i]['het_count'] == 1
            assert parsed_json['result'][i]['genes'][0] == 'TTLL5'

     

if __name__ == '__main__':
    unittest.main()

