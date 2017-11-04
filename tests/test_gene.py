
import unittest
import runserver
import sys
import load_data
import helper

class GenePageTestCase(unittest.TestCase):

    def setUp(self):
        load_data.load_data()
        runserver.app.config['TESTING'] = True
        runserver.app.config['DB_NAME'] = 'test_uclex'
        runserver.app.config['DB_NAME_HPO'] = 'test_hpo'
        runserver.app.config['DB_NAME_PATIENTS'] = 'test_patients'
        runserver.app.config['DB_NAME_USERS'] = 'test_users'
        self.app = runserver.app.test_client()
        helper.create_neo4j_demo_user()
        helper.login(self.app)

    def tearDown(self):
        self.app.get('/logout', follow_redirects=True)
    
    def gene_page(self, geneName):
        return self.app.get('/gene/'+geneName, follow_redirects=True)

    def test_gene_page(self):
        page = self.gene_page('TTLL5')
        assert page.status_code == 200
        assert 'TTLL5' in page.data 
        assert 'ENSG00000119685' in page.data
        assert 'Macular dystrophy' in page.data 
        assert 'Abnormality of the macula' in page.data
        assert 'Autosomal recessive inheritance' in page.data
        assert 'Mode of inheritance' in page.data
        assert 'Visual impairment' in page.data
        assert 'Abnormality of vision' in page.data
        assert 'Abnormal eye physiology' in page.data
        assert 'Retinal dystrophy' in page.data
        assert 'Abnormality of the retina' in page.data
        assert 'Abnormality of the fundus' in page.data
        assert 'Abnormality of the posterior segment of the globe' in page.data
        assert 'Abnormality of the globe' in page.data
        assert 'Abnormal eye morphology' in page.data
        assert 'Abnormality of the eye' in page.data
        assert 'Phenotypic abnormality' in page.data
        assert 'All' in page.data


if __name__ == '__main__':
    unittest.main()
