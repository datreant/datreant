"""
Tests for File classes.

"""

from unittest import TestCase
import MDSynthesis.Core.Files as Files
import yaml
import os

def suite():
    """Gather all tests in this module into a test suite.

    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestFile))

    return test_suite

class TestFile(TestCase):
    """Test class for the File class.

    """
    def setUp(self):
        """Create test File object. 
            
        """

        # data structure used for test Files
        self.datastruct = {'name': 'TestFile',
                           'categories': {'red': 1, 'blue': 2},
                           'tags': ['little', 'big'],
                           'number': 256}
        
        f = open('testsite/TestFile.txt', 'w')
        yaml.dump(self.datastruct, f)
        f.close()
        
        self.logger = logging.getLogger('TestFile')
        logfile = 'testsite/TestFile.log'
        fh = logging.FileHandler(logfile)
        self.logger.addHandler(fh)

    def tearDown(self):
        """Destroy File object references, and their associated files.

        """
        testfile = self.testfile.filename

        del self.testfile
        
        if os.path.exists(testfile)
            os.remove(testfile)

    def test_shlock(self):
        """Test the shared locking mechanism of the File class.

        """
        self.assertFalse(os.path.exists(self.testfile.lockname))
        self.testfile.lock()
        self.assertTrue(os.path.exists(self.testfile.lockname))

    def test_unlock(self):
        """Test the unlock mechanism of the File class.

        """
        

class TestContainerFile(TestCase):
    """Test class for ContainerFile.

    """
    def setUp(self):
        """Create test ContainerFile.

        """
        self.test = Files.ContainerFile('testcont.hdf', None, 'Sim', name='test')

    def test_create(self):
        # create a variety of files; check that results match input
        pass

    def test_add_tags(self):
        """Test out the adding of tags.

        """
        test.add_tags('mark', 932, 'mark', 'luke', '932', 8400.321)

    def test_del_tags(self):
        """Test out the deletion of tags.

        """
        test.add_tags('mark', 932, 'mark', 'luke', '932', 8400.321)
        test.del_tags('mark', 932, '932')

        # test all keyword

    def test_add_categories(self):
        """Test out adding of categories.

        """
        test.add_categories(sudo='bark', masc=29, jesus=61, xenu="lol")

    def test_del_categories(self):
        """Test out the deletion of categories.

        """
        test.add_categories(sudo='bark', masc=29, jesus=61, xenu="lol")
        test.del_categories('jesus', 'xenu', 'masc', 'xenu', all=True)

        # test all keyword



if __name__ == '__main__':
    unittest.main()

