
import unittest
import sys
import load_data
sys.path.append('offline_analysis/commons')
from phenopolis_utils import *

class utilsTestCase(unittest.TestCase):
    def setUp(self):
        self.db_config = {
            'hpo_db':'test_hpo',
            'phenopolis_db':'test_uclex',
            'patient_db':'test_patients',
            'pubmedbatch_db':'test_pubmedbatch',
        }
        load_data.load_data()
        self.db = get_mongo_collections(self.db_config)

    def tearDown(self):
        pass

    def test_mongodb(self):
        # hpo
        case = self.db['hpo_db'].hpo.find_one({'id':'HP:0000001'})
        self.assertEqual(case['name'][0],'All')

    def test_get_chrom_genes(self):
        case = get_chrom_genes([1],self.db['phenopolis_db'])
        assert 'ENSG00000196944' in case
        case = get_chrom_genes([2,1],self.db['phenopolis_db'])
        assert 'ENSG00000196944' in case

    def test_get_candidate_genes(self):
        result = get_candidate_genes(self.db['patient_db'],fields=['hpo'])
        # HP:0000548 in at least one of the returned TTLL5 result
        case = [i2 for i1 in result['ENSG00000156171']['data'] for i2 in i1['hpo'] if i2['observed'] == 'yes' and i2['id'] == 'HP:0000548']
        self.assertTrue(case)

if __name__ == '__main__':
    unittest.main()
