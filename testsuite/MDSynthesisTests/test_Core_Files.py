"""
Tests for File classes.

"""

import unittest

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
        """Create test File object.

        """

if __name__ == '__main__':
    unittest.main()

