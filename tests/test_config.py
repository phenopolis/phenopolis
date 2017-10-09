
import unittest
from config import config


class ConfigTestCase(unittest.TestCase):

    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_import_flag(self):
        assert config.IMPORT_PYSAM_PRIMER3 == True # This needs to be True before committing to the repo.    


if __name__ == '__main__':
    unittest.main()
