"""
File backends for storing pandas objects.

"""

import pandas as pd
import numpy as np


from datreant.backends.core import File


pddatafile = 'pdData.h5'


class pdDataFile(File):
    """Interface to pandas object data files.

    Data is stored as pandas data structures (Series, DataFrame, Panel) in
    the HDF5 format. This class gives the needed components for storing
    and retrieving stored data. It uses pandas' HDFStore object as its
    backend.

    """
    def _open_file_r(self):
        return pd.HDFStore(self.filename, 'r')

    def _open_file_w(self):
        return pd.HDFStore(self.filename, 'a')

    @File._write
    def add_data(self, key, data):
        """Add a pandas data object (Series, DataFrame, Panel) to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                name given to the data; used as the index for retrieving
                the data later
            *data*
                the data object to store; should be either a Series, DataFrame,
                or Panel
        """
        # index all columns if possible
        try:
            # FIXME: band-aid heuristic to catch a known corner case that
            # HDFStore doesn't catch; see ``Issue 20``
            if (isinstance(data, pd.DataFrame) and
                    data.columns.dtype == np.dtype('int64')):
                raise AttributeError

            self.handle.put(
                key, data, format='table', data_columns=True, complevel=5,
                complib='blosc')
        except AttributeError:
            self.handle.put(
                key, data, format='table', complevel=5, complib='blosc')

    @File._write
    def append_data(self, key, data):
        """Append rows to an existing pandas data object stored in the data file.

        Note that column names of new data must match those of the existing
        data. Columns cannot be appended due to the technical details of the
        HDF5 standard. To add new columns, store as a new dataset.

        :Arguments:
            *key*
                name of existing data object to append to
            *data*
                the data object whose rows are to be appended to the existing
                stored data; must have same columns (with names) as existing
                data

        """
        try:
            self.handle.append(
                key, data, data_columns=True, complevel=5, complib='blosc')
        except AttributeError:
            self.handle.append(key, data, complevel=5, complib='blosc')

    @File._read
    def get_data(self, key, **kwargs):
        """Retrieve pandas object stored in file, optionally based on where criteria.

        :Arguments:
            *key*
                name of data to retrieve

        :Keywords:
            *where*
                conditions for what rows/columns to return
            *start*
                row number to start selection
            *stop*
                row number to stop selection
            *columns*
                list of columns to return; all columns returned by default
            *iterator*
                if True, return an iterator [``False``]
            *chunksize*
                number of rows to include in iteration; implies
                ``iterator=True``

        :Returns:
            *data*
                the selected data
        """
        return self.handle.select(key, **kwargs)

    @File._write
    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        :Keywords:
            *where*
                conditions for what rows/columns to remove
            *start*
                row number to start selection
            *stop*
                row number to stop selection

        """
        self.handle.remove(key, **kwargs)

    # TODO: remove this; since we only place one datastructure in an HDF5 file,
    # we don't need it
    @File._read
    def list_data(self):
        """List names of all stored datasets.

        Although the true names start with '\' indicating the root of the
        HDF5 data tree, we wish to abstract this away. We remove the leading
        '\' from the output. This shouldn't cause any problems since the
        leading '\' can be neglected when referring to stored objects by name
        using all of pandas.HDFStore's methods anyway.

        """
        keys = self.handle.keys()
        return [i.lstrip('/') for i in keys]
