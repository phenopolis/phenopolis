
import unittest
from config import config


class ConfigTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_import_flag(self):
        assert config.IMPORT_PYSAM_PRIMER3 == True
   
    def test_phenotips_flag(self):
        assert config.LOCAL_WITH_PHENOTIPS == False     

if __name__ == '__main__':
    unittest.main()
