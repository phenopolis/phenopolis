
import unittest
import runserver

class GenePageTestCase(unittest.TestCase):

    def setUp(self):
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()

    def tearDown(self):
        pass
    
    # import this TODO LMTW
    def login(self, username, password):
        return self.app.post('/login', data=dict(
            username=username,
            password=password
        ), follow_redirects=True)

    def gene_page(self, geneName):
        return self.app.get('/gene/'+geneName, follow_redirects=True)

    def test_gene_page(self):
        rv = self.login('demo', 'demo123')
        rv = self.gene_page('TTLL5')
        assert 'TTLL5' in rv.data

if __name__ == '__main__':
    unittest.main()
