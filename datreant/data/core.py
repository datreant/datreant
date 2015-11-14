"""
Abstract interface components for reading and writing datasets.

"""

import os

import numpy as np
import pandas as pd

from datreant.data import pydata
from datreant.data import npdata
from datreant.data import pddata


class DataFile(object):
    """Interface to data files.

    This is an abstraction layer to the pddata.pdDataFile, npdata.npDataFile,
    and pydata.pyDataFile objects. This can be used by higher level objects
    without worrying about whether to use pandas storers, numpy storers, or
    pickle.

    """

    def __init__(self, datadir, logger=None, datafiletype=None, **kwargs):
        """Initialize data interface.

        :Arguments:
           *datadir*
              path to data directory
           *logger*
              Treant's logger instance
           *datafiletype*
              If known, either pddata.pddatafile or npdata.npdatafile

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
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile),
                logger=self.logger)
        elif isinstance(data, (pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D)):
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile),
                logger=self.logger)
        else:
            self.datafile = pydata.pyDataFile(
                os.path.join(self.datadir, pydata.pydatafile),
                logger=self.logger)

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
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile),
                logger=self.logger)
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
        if self.datafiletype == npdata.npdatafile:
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile),
                logger=self.logger)
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddata.pddatafile:
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile),
                logger=self.logger)
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydata.pydatafile:
            self.datafile = pydata.pyDataFile(
                os.path.join(self.datadir, pydata.pydatafile),
                logger=self.logger)
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
        if self.datafiletype == npdata.npdatafile:
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile),
                logger=self.logger)
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddata.pddatafile:
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile),
                logger=self.logger)
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydata.pydatafile:
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
        if self.datafiletype == npdata.npdatafile:
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile),
                logger=self.logger)
            out = self.datafile.list_data(key, data, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddata.pddatafile:
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile),
                logger=self.logger)
            out = self.datafile.list_data(key, data, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydata.pydatafile:
            self.logger.info('No substructure to serialized python object')
            out = None
        else:
            self.logger.info('Cannot return data without knowing datatype.')
            out = None

        return out
