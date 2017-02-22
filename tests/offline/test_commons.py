
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

    def test_gene_names_to_ids(self):
        case = gene_names_to_ids(self.db['phenopolis_db'],['LARGE','ERO1L'])
        self.assertEqual(case['LARGE']['id'],'ENSG00000133424')
        self.assertEqual(case['LARGE']['symbol'],'LARGE1')
        self.assertEqual(case['ERO1L']['id'],'ENSG00000197930')
        self.assertEqual(case['ERO1L']['symbol'],'ERO1A')

    def test_get_candidate_genes(self):
        result = get_candidate_genes(self.db,fields=['hpo'])
        case = result['ENSG00000156171']
        self.assertEqual(case['symbol'],'DRAM2')
        result = get_candidate_genes(self.db,genes=['ABCA4','DRAM2'],fields=['hpo'])
        case = result['ENSG00000156171']
        self.assertEqual(case['symbol'],'DRAM2')
        

if __name__ == '__main__':
    unittest.main()
