"""
Tests for File classes.

"""

import unittest
import MDSynthesis.Core.Files
import yaml
import pickle as pkl
import json

def suite():
    """Gather all tests in this module into a test suite.

    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestFile))

    return test_suite

class TestFile(unittest.TestCase):
    """Test class for the File class.

    """
    def setUp(self):
        """Create test File objects. 
            
        We wish to test using multiple readers and writers at the same time. 

        """

        # data structure used for test Files
        self.datastruct = {'name': 'TestFile',
                           'categories': {'red': 1, 'blue': 2},
                           'tags': ['little', 'big'],
                           'number': 256}
        
        self.testfile_yaml = MDSynthesis.Core.Files.File(filename='testsite/TestFile.yaml',
                                reader=yaml.load, writer=yaml.dump, datastruct=self.datastruct}
        
        self.logger_yaml = logging.getLogger('TestFileYaml')
        logfile = 'testsite/TestFileYaml.log'
        fh = logging.FileHandler(logfile)
        self.logger_yaml.addHandler(fh)


        self.testfile_pkl = MDSynthesis.Core.Files.File(filename='testsite/TestFile.pkl',
                                reader=pkl.load, writer=pkl.dump, datastruct=self.datastruct}

        self.logger_pkl = logging.getLogger('TestFilePkl')
        logfile = 'testsite/TestFilePkl.log'
        fh = logging.FileHandler(logfile)
        self.logger_pkl.addHandler(fh)


        self.testfile_json = MDSynthesis.Core.Files.File(filename='testsite/TestFile.json',
                                reader=json.load, writer=json.dump, datastruct=self.datastruct}

        self.logger_json = logging.getLogger('TestFileJson')
        logfile = 'testsite/TestFileJson.log'
        fh = logging.FileHandler(logfile)
        self.logger_json.addHandler(fh)


    def test_lock(self):
        """Test the basic locking mechanism of the File class.

        """
        pass
        

if __name__ == '__main__':
    unittest.main()

