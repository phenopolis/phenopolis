
import unittest
import sys
sys.path.append('offline_analysis/pubmedScore')
import pubmedScore
import time

class pubmedScoreTestCase(unittest.TestCase):
    def setUp(self):
        keywords = "retina,retinal,retinitis,blindness,macula,macular,stargardt,pigmentosa"
        self.keywords = keywords.split(',')
        self.now = time.mktime(time.localtime())

    def tearDown(self):
        pass

    def test_find_item(self):
        obj = {'foo':1}
        case = pubmedScore.find_item(obj,'foo')
        self.assertEqual(case,1)
        case = pubmedScore.find_item(obj,'fo')
        self.assertEqual(case,None)
        
        obj = {'bar':{'foo':'haha'},'baz':0}
        case = pubmedScore.find_item(obj,'foo')
        self.assertEqual(case,'haha')

    def test_score(self):
        # case might change
        case = pubmedScore.pubmed_query('wahahahaha0>o', self.keywords)
        self.assertEqual(case['total_score'],0)
        case = pubmedScore.pubmed_query('ARL2BP', self.keywords)
        c = [i for i in case['results'] if i['id'] == '23849777']
        self.assertTrue(c)

    def test_pubmed(self):
        case = pubmedScore.pubmed('wahahahaha0>o', self.keywords, self.now, test=True)
        self.assertEqual(case['total_score'],0)
        case = pubmedScore.pubmed('TTLL5', self.keywords, self.now, test=True)
        c = [i for i in case['results'] if i['id'] == '27162334']
        self.assertTrue(c)

if __name__ == '__main__':
    unittest.main()
