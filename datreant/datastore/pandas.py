"""
Interface classes for state files and data files.

"""

import os
import sys
import fcntl
import pickle
import logging
import warnings
from functools import wraps

import tables
import yaml
import h5py
import pandas as pd
import numpy as np

import datreant

# pandas Datafile
pddatafile = "pdData.h5"

# numpy Datafile
npdatafile = "npData.h5"

# catchall DataFile
pydatafile = "pyData.pkl"

class DataFile(object):
    """Interface to data files.

    This is an abstraction layer to the pdDataFile, npDataFile, and pyDataFile
    objects. This can be used by higher level objects without worrying about
    whether to use pandas storers, numpy storers, or pickle.

    """

    def __init__(self, datadir, logger=None, datafiletype=None, **kwargs):
        """Initialize data interface.

        :Arguments:
           *datadir*
              path to data directory
           *logger*
              Treant's logger instance
           *datafiletype*
              If known, either pddatafile or npdatafile

        """
        self.datadir = datadir
        self.datafile = None

        # if given, can get data
        self.datafiletype = datafiletype

        self.logger = logger

    def add_data(self, key, data):
        """Add a pandas data object (Series, DataFrame, Panel), numpy array,
        or pickleable python object to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                name given to the data; used as the index for retrieving
                the data later
            *data*
                the data object to store; should be either a pandas Series,
                DataFrame, Panel, or a numpy array
        """
        if isinstance(data, np.ndarray):
            self.datafile = npDataFile(
                os.path.join(self.datadir, npdatafile), logger=self.logger)
        elif isinstance(data, (pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D)):
            self.datafile = pdDataFile(
                os.path.join(self.datadir, pddatafile), logger=self.logger)
        else:
            self.datafile = pyDataFile(
                os.path.join(self.datadir, pydatafile), logger=self.logger)

        self.datafile.add_data(key, data)

        # dereference
        self.datafile = None

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
        # TODO: add exceptions where appending isn't possible
        if isinstance(data, np.ndarray):
            self.logger.info('Cannot append numpy arrays.')
        elif isinstance(data, (pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D)):
            self.datafile = pdDataFile(
                os.path.join(self.datadir, pddatafile), logger=self.logger)
            self.datafile.append_data(key, data)
        else:
            self.logger.info('Cannot append python object.')

            # dereference
            self.datafile = None

    def get_data(self, key, **kwargs):
        """Retrieve data object stored in file.

        :Arguments:
            *key*
                name of data to retrieve

        :Keywords:
            *where*
                for pandas objects, conditions for what rows/columns to return
            *start*
                for pandas objects, row number to start selection
            *stop*
                for pandas objects, row number to stop selection
            *columns*
                for pandas objects, list of columns to return; all columns
                returned by default
            *iterator*
                for pandas objects, if True, return an iterator [``False``]
            *chunksize*
                for pandas objects, number of rows to include in iteration;
                implies ``iterator=True``

        :Returns:
            *data*
                the selected data
        """
        if self.datafiletype == npdatafile:
            self.datafile = npDataFile(
                os.path.join(self.datadir, npdatafile), logger=self.logger)
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddatafile:
            self.datafile = pdDataFile(
                os.path.join(self.datadir, pddatafile), logger=self.logger)
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydatafile:
            self.datafile = pyDataFile(
                os.path.join(self.datadir, pydatafile), logger=self.logger)
            out = self.datafile.get_data(key)
            self.datafile = None
        else:
            # TODO: add exception here
            self.logger.info('Cannot return data without knowing datatype.')
            out = None

        return out

    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        :Keywords:
            *where*
                for pandas objects, conditions for what rows/columns to remove
            *start*
                for pandas objecs, row number to start selection
            *stop*
                for pandas objects, row number to stop selection

        """
        if self.datafiletype == npdatafile:
            self.datafile = npDataFile(
                os.path.join(self.datadir, npdatafile), logger=self.logger)
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddatafile:
            self.datafile = pdDataFile(
                os.path.join(self.datadir, pddatafile), logger=self.logger)
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydatafile:
            pass
        else:
            # TODO: add exception here
            self.logger.info('Cannot return data without knowing datatype.')
            out = None

    # TODO: remove this; since we only place one datastructure in an HDF5 file,
    # we don't need it
    def list_data(self):
        """List names of all stored datasets.

        Although the true names start with '\' indicating the root of the
        HDF5 data tree, we wish to abstract this away. We remove the leading
        '\' from the output. This shouldn't cause any problems since the
        leading '\' can be neglected when referring to stored objects by name
        using all of pandas.HDFStore and h5py.File methods anyway.

        """
        if self.datafiletype == npdatafile:
            self.datafile = npDataFile(
                os.path.join(self.datadir, npdatafile), logger=self.logger)
            out = self.datafile.list_data(key, data, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddatafile:
            self.datafile = pdDataFile(
                os.path.join(self.datadir, pddatafile), logger=self.logger)
            out = self.datafile.list_data(key, data, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydatafile:
            self.logger.info('No substructure to serialized python object')
            out = None
        else:
            self.logger.info('Cannot return data without knowing datatype.')
            out = None

        return out


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


class npDataFile(File):
    """Interface to numpy object data files.

    Data is stored as numpy arrays in the HDF5 format. This class gives the
    needed components for storing and retrieving stored data. It uses h5py as
    its backend.

    """
    def _open_file_r(self):
        return h5py.File(self.filename, 'r')

    def _open_file_w(self):
        return h5py.File(self.filename, 'a')

    @File._write
    def add_data(self, key, data):
        """Add a numpy array to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                name given to the data; used as the index for retrieving
                the data later
            *data*
                the numpy array to store
        """
        try:
            self.handle.create_dataset(key, data=data)
        except RuntimeError:
            del self.handle[key]
            self.handle.create_dataset(key, data=data)

    @File._read
    def get_data(self, key, **kwargs):
        """Retrieve numpy array stored in file.

        :Arguments:
            *key*
                name of data to retrieve

        :Returns:
            *data*
                the selected data
        """
        return self.handle[key].value

    @File._write
    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        """
        del self.handle[key]

    @File._read
    def list_data(self):
        """List names of all stored datasets.

        Although the true names start with '\' indicating the root of the
        HDF5 data tree, we wish to abstract this away. We remove the leading
        '\' from the output. This shouldn't cause any problems since the
        leading '\' can be neglected when referring to stored objects by name
        using all of h5py.File's methods anyway.

        """
        keys = self.handle.keys()
        return keys


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
