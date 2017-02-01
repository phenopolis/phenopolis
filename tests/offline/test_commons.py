
import unittest
import sys
sys.path.append('offline_analysis/commons')
sys.path.append('tests')
import load_data
from phenopolis_utils import *

class utilsTestCase(unittest.TestCase):
    def setUp(self):
        self.db_config = {
            'hpo_db':'test_hpo',
            'phenopolis_db':'test_uclex',
            'patient_db':'test_patients',
        }

    def test_mongodb(self):
        # hpo
        db = get_mongo_collections(self.db_config)
        case = db['hpo_db'].hpo.find_one({'id':'HP:0000001'})
        self.assertEqual(case['name'][0],'All')

if __name__ == '__main__':
    unittest.main()
