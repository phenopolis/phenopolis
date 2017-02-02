
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

if __name__ == '__main__':
    unittest.main()
