"""
File backends for storing general python objects.

"""

import pickle

from datreant.backends.core import File


pydatafile = 'pyData.h5'


class pyDataFile(File):
    """Interface to python object data files.

    Arbitrary python objects are stored as pickled objects on disk. This class
    gives the needed components for storing and retrieving stored data in the
    same basic way as for pandas and numpy objects. It uses pickle files for
    serialization.

    """
    def _open_file_r(self):
        return open(self.filename, 'rb')

    def _open_file_w(self):
        return open(self.filename, 'wb+')

    @File._write
    def add_data(self, key, data):
        """Add a numpy array to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                not used, but needed to give consistent interface
            *data*
                the numpy array to store
        """
        pickle.dump(data, self.handle)

    @File._read
    def get_data(self, key, **kwargs):
        """Retrieve numpy array stored in file.

        :Arguments:
            *key*
                not used, but needed to give consistent interface

        :Returns:
            *data*
                the selected data
        """
        return pickle.load(self.handle)
