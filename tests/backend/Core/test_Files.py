"""
Tests for File classes.

"""

import unittest
from unittest import TestCase
import MDSynthesis.Core.Files as Files
import os

def suite():
    """Gather all tests in this module into a test suite.

    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestFile))
    test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContainerFile))

    return test_suite

class TestFile(TestCase):
    """Test class for the File class.

    """
    def setUp(self):
        """Create test File object. 
            
        """
        pass

    def tearDown(self):
        """Destroy File object references, and their associated files.

        """
        pass

    def test_shlock(self):
        """Test the shared locking mechanism of the File class.

        """
        pass

    def test_unlock(self):
        """Test the unlock mechanism of the File class.

        """
        pass
        
class TestContainerFile(TestCase):
    """Test class for ContainerFile.

    """
    def setUp(self):
        """Create test ContainerFile.

        """
        self.test = Files.ContainerFile('testcont.hdf', None, 'Sim', name='test')

    def tearDown(self):
        """Delete test ContainerFile.

        """
        os.remove('testcont.hdf')

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
        inputs = {'sudo': 'bark', 'masc': 29, 'jesus': 61, 'xenu': 'lol'}

        test.add_categories(**inputs)
        test.del_categories('jesus', 'xenu', 'masc', 'xenu')

        test._open_r()
        table = test.handle.root.categories
        result = { x['category']: x['value'] for x in table.iterrows() }
        test._close()

        self.assertEqual(result, {'sudo': 'bark'})

    def test_update_name(self):
        """Test rename of Container.

        """
        newname = 'blark'
        test.update_name(newname)

        test._open_r()
        table = test.handle.root.meta
        result = table.cols.name[0]
        test._close()

        self.assertEqual(newname, result)

if __name__ == '__main__':
    unittest.main()

