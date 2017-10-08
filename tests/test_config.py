
import unittest
from config import config


class ConfigTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_import_flag(self):
        assert config.IMPORT_PYSAM_PRIMER3 == True # This needs to be True before committing to the repo.    

    # This needs to be False before committing to the repo because the code 
    # that uses py2neo is not set up to work with Travis.    
    def test_py2neo_flag(self):
        assert config.USE_PY2NEO == False 


if __name__ == '__main__':
    unittest.main()
