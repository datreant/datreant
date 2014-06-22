"""
Tests for File classes.

"""

import unittest
import MDSynthesis.Core.Files
import yaml
import os

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
        For now we just test with yaml.

        """

        # data structure used for test Files
        self.datastruct = {'name': 'TestFile',
                           'categories': {'red': 1, 'blue': 2},
                           'tags': ['little', 'big'],
                           'number': 256}
        
        self.testfile = MDSynthesis.Core.Files.File(filename='testsite/TestFile.yaml',
                             reader=yaml.load, writer=yaml.dump, datastruct=self.datastruct}
        
        self.logger = logging.getLogger('TestFile')
        logfile = 'testsite/TestFile.log'
        fh = logging.FileHandler(logfile)
        self.logger.addHandler(fh)

    def tearDown(self):
        """

        """


    def test_lock(self):
        """Test the basic locking mechanism of the File class.

        """
        self.testfile.lock()
        self.assertTrue(os.path.exists(self.testfile.lockname))


if __name__ == '__main__':
    unittest.main()

