
import unittest
import runserver
import load_data

class GenePageTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        runserver.app.config['DB_NAME'] = 'test_uclex'
        runserver.app.config['DB_NAME_HPO'] = 'test_hpo'
        runserver.app.config['DB_NAME_PATIENTS'] = 'test_patients'
        self.app = runserver.app.test_client()

    def tearDown(self):
        pass
    
    # TODO LMTW - Make this a common helper function 
    def login(self, username, password):
        return self.app.post('/login', data=dict(
            username=username,
            password=password
        ), follow_redirects=True)

    def gene_page(self, geneName):
        return self.app.get('/gene/'+geneName, follow_redirects=True)

    def test_gene_page(self):
        load_data.load_data()
        rv = self.login('demo', 'demo123')
        rv = self.gene_page('TTLL5')
        assert 'TTLL5' in rv.data 
        assert 'Macular dystrophy' in rv.data 
        assert 'Abnormality of the retina' in rv.data
        assert 'ENSG00000119685' in rv.data

if __name__ == '__main__':
    unittest.main()
