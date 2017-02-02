
import unittest
import sys
sys.path.append('offline_analysis/pubmedScore')
import pubmedScore
import time

KEYWORDS = "retina,retinal,retinitis,blindness,macula,macular,stargardt,pigmentosa"
KEYWORDS = KEYWORDS.split(',')
NOW = time.mktime(time.localtime())

class pubmedScoreTestCase(unittest.TestCase):
    def test_score(self):
        # case might change
        case = pubmedScore.pubmed_query('wahahahaha0>o', KEYWORDS)
        self.assertEqual(case['total_score'],0)
        case = pubmedScore.pubmed_query('ARL2BP', KEYWORDS)
        c = [i for i in case['results'] if i['id'] == '23849777']
        self.assertTrue(c)
    def test_pubmed(self):
        case = pubmedScore.pubmed('wahahahaha0>o', KEYWORDS, NOW, test=True)
        self.assertEqual(case['total_score'],0)
        case = pubmedScore.pubmed('TTLL5', KEYWORDS, NOW, test=True)
        c = [i for i in case['results'] if i['id'] == '27162334']
        self.assertTrue(c)

if __name__ == '__main__':
    unittest.main()
