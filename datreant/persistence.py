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
import h5py
import pandas as pd
import numpy as np

# extension used for state files
statefile_ext = 'h5'

# pandas Datafile
pddatafile = "pdData.h5"

# numpy Datafile
npdatafile = "npData.h5"

# catchall DataFile
pydatafile = "pyData.pkl"

# number of characters required for uuids
uuidlength = 36

# max length in characters for all paths
pathlength = 511

# max character length of strings used for handles, tags, categories
namelength = 55


# TODO: add careful checks that file is actually a state file
# TODO: must be extendable in derived packages to give Treant subclasses
def treantfile(filename, logger=None, **kwargs):
    """Generate or regenerate the appropriate treant file instance from
    filename.

    :Arguments:
        *filename*
            path to state file (existing or to be created), including the
            filename
        *logger*
            logger instance to pass to treant file instance

    **kwargs passed to treant file ``__init__()`` method

    """
    basename = os.path.basename(filename)
    if 'Treant' in basename:
        statefileclass = TreantFile
    elif 'Group' in basename:
        statefileclass = GroupFile

    return statefileclass(filename, logger=logger, **kwargs)


class File(object):
    """File object base class. Implements file locking and reloading methods.

    """

    def __init__(self, filename, logger=None, **kwargs):
        """Create File instance for interacting with file on disk.

        All files in MDSynthesis should be accessible by high-level methods
        without having to worry about simultaneous reading and writing by other
        processes. The File object includes methods and infrastructure for
        ensuring shared and exclusive locks are consistently applied before
        reads and writes, respectively. It handles any other low-level tasks
        for maintaining file integrity.

        :Arguments:
            *filename*
                name of file on disk object corresponds to
            *logger*
                logger to send warnings and errors to

        """
        # filter NaturalNameWarnings from pytables, when they arrive
        warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

        self.filename = os.path.abspath(filename)
        self.handle = None

        self._start_logger(logger)

    def _start_logger(self, logger):
        """Start up the logger.

        """
        # delete current logger
        try:
            del self.logger
        except AttributeError:
            pass

        # log to standard out if no logger given
        if not logger:
            self.logger = logging.getLogger(
                '{}'.format(self.__class__.__name__))
            self.logger.setLevel(logging.INFO)

            if not any([isinstance(x, logging.StreamHandler)
                        for x in self.logger.handlers]):
                ch = logging.StreamHandler(sys.stdout)
                cf = logging.Formatter(
                        '%(name)-12s: %(levelname)-8s %(message)s')
                ch.setFormatter(cf)
                self.logger.addHandler(ch)
        else:
            self.logger = logger

    def _shlock(self, f):
        """Get shared lock on file.

        Using fcntl.lockf, a shared lock on the file is obtained. If an
        exclusive lock is already held on the file by another process,
        then the method waits until it can obtain the lock.

        :Arguments:
            *f*
                file handle

        :Returns:
            *success*
                True if shared lock successfully obtained
        """
        fcntl.lockf(f, fcntl.LOCK_SH)

        return True

    def _exlock(self, f):
        """Get exclusive lock on file.

        Using fcntl.lockf, an exclusive lock on the file is obtained. If a
        shared or exclusive lock is already held on the file by another
        process, then the method waits until it can obtain the lock.

        :Arguments:
            *f*
                file handle

        :Returns:
            *success*
                True if exclusive lock successfully obtained
        """
        fcntl.lockf(f, fcntl.LOCK_EX)

        return True

    def _unlock(self, f):
        """Remove exclusive or shared lock on file.

        WARNING: It is very rare that this is necessary, since a file must be
        unlocked before it is closed. Furthermore, locks disappear when a file
        is closed anyway.  This method will remain here for now, but may be
        removed in the future if not needed (likely).

        :Arguments:
            *f*
                file handle

        :Returns:
            *success*
                True if lock removed
        """
        fcntl.lockf(f, fcntl.LOCK_UN)

        return True

    def _check_existence(self):
        """Check for existence of file.

        """
        return os.path.exists(self.filename)

    @staticmethod
    def _read_state(func):
        """Decorator for opening state file for reading and applying shared
        lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.handle.isopen:
                out = func(self, *args, **kwargs)
            else:
                self.handle = tables.open_file(self.filename, 'r')
                self._shlock(self.handle)
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
            return out

        return inner

    @staticmethod
    def _write_state(func):
        """Decorator for opening state file for writing and applying exclusive lock.

        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method
        is executed. It also ensures that the lock is removed and the file
        closed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            self.handle = tables.open_file(self.filename, 'a')
            self._exlock(self.handle)
            try:
                out = func(self, *args, **kwargs)
            finally:
                self.handle.close()
            return out

        return inner


class TreantFile(File):
    """Treant file object; syncronized access to Treant data.

    """
    class _Meta(tables.IsDescription):
        """Table definition for metadata.

        All strings limited to hardcoded size for now.

        """
        # version of MDSynthesis object file data corresponds to
        # allows future-proofing of old objects so that formats of new releases
        # can be automatically built from old ones
        version = tables.StringCol(15)

    class _Coordinator(tables.IsDescription):
        """Table definition for coordinator info.

        This information is kept separate from other metadata to allow the
        Coordinator to simply stack tables to populate its database. It doesn't
        need entries that store its own path.

        Path length fixed size for now.
        """
        # absolute path of coordinator
        abspath = tables.StringCol(pathlength)

    class _Tags(tables.IsDescription):
        """Table definition for tags.

        """
        tag = tables.StringCol(namelength)

    class _Categories(tables.IsDescription):
        """Table definition for categories.

        """
        category = tables.StringCol(namelength)
        value = tables.StringCol(namelength)

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Treant state file.

        This is the base class for all Treant state files. It generates data
        structure elements common to all Treants. It also implements
        low-level I/O functionality.

        :Arguments:
            *filename*
                path to file
            *logger*
                Treant's logger instance

        :Keywords:
            *treanttype*
                Treant type
            *name*
                user-given name of Treant object
            *coordinator*
                directory in which coordinator state file can be found [None]
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search
            *version*
                version of MDSynthesis file was generated with

        .. Note:: kwargs passed to :meth:`create`

        """
        super(TreantFile, self).__init__(filename, logger=logger)

        # if file does not exist, it is created
        if not self._check_existence():
            self.create(**kwargs)
        else:
            # open file for the first time to initialize handle
            self.handle = tables.open_file(self.filename, 'r')
            self.handle.close()

    def create(self, **kwargs):
        """Build state file and common data structure elements.

        :Keywords:
            *treanttype*
                Treant type
            *name*
                user-given name of Treant object
            *coordinator*
                directory in which coordinator state file can be found [None]
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search
            *version*
                version of MDSynthesis file was generated with
        """
        import datreant

        treanttype = kwargs.pop('treanttype', None)

        # metadata table
        self.update_version(kwargs.pop('version', datreant.__version__))

        # coordinator table
        self.update_coordinator(kwargs.pop('coordinator', None))

        # tags table
        tags = kwargs.pop('tags', list())
        self.add_tags(*tags)

        # categories table
        categories = kwargs.pop('categories', dict())
        self.add_categories(**categories)

    def get_location(self):
        """Get Treant basedir.

        :Returns:
            *location*
                absolute path to Treant basedir

        """
        return os.path.dirname(self.filename)

    @File._read_state
    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.version[0]

    @File._write_state
    def update_version(self, version):
        """Update version of Treant.

        :Arugments:
            *version*
                new version of Treant

        """
        try:
            table = self.handle.get_node('/', 'meta')
            table.cols.version[0] = version
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'meta', self._Meta, 'metadata')
            table.row['version'] = version
            table.row.append()

    @File._read_state
    def get_coordinator(self):
        """Get absolute path to Coordinator.

        :Returns:
            *coordinator*
                absolute path to Coordinator directory

        """
        table = self.handle.get_node('/', 'coordinator')
        out = table.cols.abspath[0]

        if out == 'None':
            out = None
        return out

    @File._write_state
    def update_coordinator(self, coordinator):
        """Update Treant location.

        :Arguments:
            *coordinator*
                absolute path to Coordinator directory
        """
        try:
            table = self.handle.get_node('/', 'coordinator')
            if coordinator:
                table.cols.abspath[0] = os.path.abspath(coordinator)
            else:
                table.cols.abspath[0] = None
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'coordinator', self._Coordinator,
                'coordinator information')
            if coordinator:
                table.row['abspath'] = os.path.abspath(coordinator)
            else:
                table.row['abspath'] = None
            table.row.append()

    @File._read_state
    def get_tags(self):
        """Get all tags as a list.

        :Returns:
            *tags*
                list of all tags
        """
        table = self.handle.get_node('/', 'tags')
        return [x['tag'] for x in table.read()]

    @File._write_state
    def add_tags(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
        one another. Sometimes preferable to categories.

        :Arguments:
            *tags*
                Tags to add. Must be convertable to strings using the str()
                builtin.

        """
        try:
            table = self.handle.get_node('/', 'tags')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'tags', self._Tags, 'tags')

        # ensure tags are unique (we don't care about order)
        tags = set([str(tag) for tag in tags])

        # remove tags already present in metadata from list
        # TODO: more efficient way to do this?
        tags_present = list()
        for row in table:
            for tag in tags:
                if (row['tag'] == tag):
                    tags_present.append(tag)

        tags = list(tags - set(tags_present))

        # add new tags
        for tag in tags:
            table.row['tag'] = tag
            table.row.append()

    @File._write_state
    def del_tags(self, *tags, **kwargs):
        """Delete tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.

        :Keywords:
            *all*
                When True, delete all tags [``False``]

        """
        table = self.handle.get_node('/', 'tags')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table('/', 'tags', self._Tags, 'tags')

        else:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])

            # get matching rows
            rowlist = list()
            for row in table:
                for tag in tags:
                    if (row['tag'] == tag):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'tags', self._Tags, 'tags')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1

    @File._read_state
    def get_categories(self):
        """Get all categories as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories
        """
        table = self.handle.get_node('/', 'categories')
        return {x['category']: x['value'] for x in table.read()}

    @File._write_state
    def add_categories(self, **categories):
        """Add any number of categories to the Treant.

        Categories are key-value pairs of strings that serve to differentiate
        Treants from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        :Keywords:
            *categories*
                Categories to add. Keyword used as key, value used as value.
                Both must be convertible to strings using the str() builtin.

        """
        try:
            table = self.handle.get_node('/', 'categories')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'categories', self._Categories, 'categories')

        table = self.handle.get_node('/', 'categories')

        # remove categories already present in metadata from dictionary
        # TODO: more efficient way to do this?
        for row in table:
            for key in categories.keys():
                if (row['category'] == key):
                    row['value'] = str(categories[key])
                    row.update()
                    # dangerous? or not since we are iterating through
                    # categories.keys() and not categories?
                    categories.pop(key)

        # add new categories
        for key in categories.keys():
            table.row['category'] = key
            table.row['value'] = str(categories[key])
            table.row.append()

    @File._write_state
    def del_categories(self, *categories, **kwargs):
        """Delete categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        :Keywords:
            *all*
                When True, delete all categories [``False``]

        """
        table = self.handle.get_node('/', 'categories')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table(
                '/', 'categories', self._Categories, 'categories')
        else:
            # remove redundant categories from given list if present
            categories = set([str(category) for category in categories])

            # get matching rows
            rowlist = list()
            for row in table:
                for category in categories:
                    if (row['category'] == category):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'categories', self._Categories, 'categories')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1

    def _open_r(self):
        """Open file with intention to write.

        Not to be used except for debugging files.

        """
        self.handle = tables.open_file(self.filename, 'r')
        self._shlock(self.handle)

    def _open_w(self):
        """Open file with intention to write.

        Not to be used except for debugging files.

        """
        self.handle = tables.open_file(self.filename, 'a')
        self._exlock(self.handle)

    def _close(self):
        """Close file.

        Not to be used except for debugging files.

        """
        self.handle.close()


class GroupFile(TreantFile):
    """Main Group state file.

    This file contains all the information needed to store the state of a
    Group object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.

    """
    # add new paths to include them in member searches
    memberpaths = ['abspath', 'relCont']

    class _Members(tables.IsDescription):

        """Table definition for the members of the Group.

        Stores for each member its treant type, uuid, and two versions of
        the path to the member treant: the absolute path (abspath) and the
        relative path from the Group object's directory (relCont). This allows
        the Group object to use some heuristically good starting points when
        trying to find missing files using a Foxhound.

        """
        # unique identifier for treant
        uuid = tables.StringCol(uuidlength)

        # treant type
        treanttype = tables.StringCol(namelength)

        abspath = tables.StringCol(pathlength)
        relCont = tables.StringCol(pathlength)

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Group state file.

        :Arguments:
           *filename*
              path to file
           *logger*
              logger to send warnings and errors to

        :Keywords:
           *name*
              user-given name of Treant object
           *coordinator*
              directory in which coordinator state file can be found [None]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search
        """
        super(GroupFile, self).__init__(filename, logger=logger, **kwargs)

    def create(self, **kwargs):
        """Build Group data structure.

        :Keywords:
           *name*
              user-given name of Group object
           *coordinator*
              directory in which Coordinator state file can be found [``None``]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search

        .. Note:: kwargs passed to :meth:`create`

        """
        super(GroupFile, self).create(treanttype='Group', **kwargs)

        self._make_membertable()

    @File._write_state
    def _make_membertable(self):
        """Make member table.

        Used only on file creation.

        """
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

    @File._write_state
    def add_member(self, uuid, treanttype, basedir):
        """Add a member to the Group.

        If the member is already present, its basedir paths will be updated
        with the given basedir.

        :Arguments:
            *uuid*
                the uuid of the new member
            *treanttype*
                the treant type of the new member
            *basedir*
                basedir of the new member in the filesystem

        """
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

        # check if uuid already present
        rownum = [row.nrow for row in table.where("uuid=='{}'".format(uuid))]
        if rownum:
            # if present, update location
            table.cols.abspath[rownum[0]] = os.path.abspath(basedir)
            table.cols.relCont[rownum[0]] = os.path.relpath(
                    basedir, self.get_location())
        else:
            table.row['uuid'] = uuid
            table.row['treanttype'] = treanttype
            table.row['abspath'] = os.path.abspath(basedir)
            table.row['relCont'] = os.path.relpath(
                    basedir, self.get_location())
            table.row.append()

    @File._write_state
    def del_member(self, *uuid, **kwargs):
        """Remove a member from the Group.

        :Arguments:
            *uuid*
                the uuid(s) of the member(s) to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        table = self.handle.get_node('/', 'members')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uid) for uid in uuid])

            # get matching rows
            # TODO: possibly faster to use table.where
            rowlist = list()
            for row in table:
                for uuid in uuids:
                    if (row['uuid'] == uuid):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'members', self._Members, 'members')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1

    @File._read_state
    def get_member(self, uuid):
        """Get all stored information on the specified member.

        Returns a dictionary whose keys are column names and values the
        corresponding values for the member.

        :Arguments:
            *uuid*
                uuid of the member to retrieve information for

        :Returns:
            *memberinfo*
                a dictionary containing all information stored for the
                specified member
        """
        table = self.handle.get_node('/', 'members')
        fields = table.dtype.names

        memberinfo = None
        for row in table.where("uuid == '{}'".format(uuid)):
            memberinfo = row.fetch_all_fields()

        if memberinfo:
            memberinfo = {x: y for x, y in zip(fields, memberinfo)}

        return memberinfo

    @File._read_state
    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                structured array giving full member data, with
                each row corresponding to a member
        """
        table = self.handle.get_node('/', 'members')
        return table.read()

    @File._read_state
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                array giving treanttype of each member, in order
        """
        table = self.handle.get_node('/', 'members')
        return table.read()['uuid']

    @File._read_state
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                array giving treanttype of each member, in order
        """
        table = self.handle.get_node('/', 'members')
        return table.read()['treanttype']

    @File._read_state
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                structured array containing all paths to member basedirs
        """
        table = self.handle.get_node('/', 'members')
        return table.read()[self.memberpaths]


# TODO: replace use of this class with a function that returns the proper data
# file
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

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize data file.

        :Arguments:
           *filename*
              path to file
           *logger*
              Treant's logger instance

        """
        super(pdDataFile, self).__init__(filename, logger=logger)

        # open file for the first time to initialize handle
        self.handle = pd.HDFStore(self.filename, 'a')
        self.handle.close()

    def _read_pddata(func):
        """Decorator for opening data file for reading and applying shared lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.handle.is_open:
                out = func(self, *args, **kwargs)
            else:
                self.handle = pd.HDFStore(self.filename, 'r')
                self._shlock(self.handle._handle)
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
            return out

        return inner

    def _write_pddata(func):
        """Decorator for opening data file for writing and applying exclusive lock.

        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method
        is executed. It also ensures that the lock is removed and the file
        closed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            self.handle = pd.HDFStore(self.filename, 'a')
            self._exlock(self.handle._handle)
            try:
                out = func(self, *args, **kwargs)
            finally:
                self.handle.close()
            return out

        return inner

    @_write_pddata
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

    @_write_pddata
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

    @_read_pddata
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

    @_write_pddata
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
    @_read_pddata
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

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize data file.

        :Arguments:
           *filename*
              path to file
           *logger*
              Treant's logger instance

        """
        super(npDataFile, self).__init__(filename, logger=logger)

        # open file for the first time to initialize handle
        self.handle = h5py.File(self.filename, 'a')
        self.handle.close()

    def _read_npdata(func):
        """Decorator for opening data file for reading and applying shared lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            try:
                self.handle.mode
                out = func(self, *args, **kwargs)
            except ValueError:
                self.handle = h5py.File(self.filename, 'r')
                self._shlock(self.handle.fid.get_vfd_handle())
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
            return out

        return inner

    def _write_npdata(func):
        """Decorator for opening data file for writing and applying exclusive lock.

        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method
        is executed. It also ensures that the lock is removed and the file
        closed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            self.handle = h5py.File(self.filename, 'a')
            self._exlock(self.handle.fid.get_vfd_handle())
            try:
                out = func(self, *args, **kwargs)
            finally:
                self.handle.close()
            return out

        return inner

    @_write_npdata
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

    @_read_npdata
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

    @_write_npdata
    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        """
        del self.handle[key]

    @_read_npdata
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

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize data file.

        :Arguments:
           *filename*
              path to file
           *logger*
              Treant's logger instance

        """
        super(pyDataFile, self).__init__(filename, logger=logger)

        # open file for the first time to initialize handle
        self.handle = open(self.filename, 'ab+')
        self.handle.close()

    def _read_pydata(func):
        """Decorator for opening data file for reading and applying shared lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            try:
                self.handle.fileno()
                out = func(self, *args, **kwargs)
            except ValueError:
                self.handle = open(self.filename, 'rb')
                self._shlock(self.handle)
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
            return out

        return inner

    def _write_pydata(func):
        """Decorator for opening data file for writing and applying exclusive lock.

        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method
        is executed. It also ensures that the lock is removed and the file
        closed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            self.handle = open(self.filename, 'wb+')
            self._exlock(self.handle)
            try:
                out = func(self, *args, **kwargs)
            finally:
                self.handle.close()
            return out

        return inner

    @_write_pydata
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

    @_read_pydata
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
