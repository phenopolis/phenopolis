
import unittest
import sys
sys.path.append('offline_analysis/commons')
from phenopolis_utils import *


class utilsTestCase(unittest.TestCase):
    def test_mongodb(self):
        # hpo
        db = get_mongo_collections()
        case = db['hpo_db'].hpo.find_one({'id':'HP:0000001'})
        self.assertEqual(case['name'][0],'All')

if __name__ == '__main__':
    unittest.main()
