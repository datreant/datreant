"""
Classes for datafile syncronization. 

"""

import Aggregators
import Workers
from uuid import uuid4
import tables
import pandas
import fcntl
import os
import sys
import logging
from functools import wraps

class File(object):
    """File object base class. Implements file locking and reloading methods.
    
    """
    def __init__(self, filename, logger=None, **kwargs):
        """Create File instance for interacting with file on disk.

        All files in MDSynthesis should be accessible by high-level
        methods without having to worry about simultaneous reading and writing by
        other processes. The File object includes methods and infrastructure
        for ensuring shared and exclusive locks are consistently applied before
        reads and writes, respectively. It handles any other low-level tasks
        for maintaining file integrity.

        :Arguments:
           *filename*
              name of file on disk object corresponds to 
           *logger*
              logger to send warnings and errors to

        """
        self.filename = os.path.abspath(filename)
        self.handle = None

        # log to standard out if no logger given
        if not logger:
            self.logger = logging.getLogger('{}'.format(self.__class__.__name__))
            self.logger.setLevel(logging.INFO)

            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self.logger.addHandler(ch)
        else:
            self.logger = logger

    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)
    
            # file handler
            if 'basedir' in self.metadata:
                logfile = os.path.join(self.metadata['basedir'], self._containerlog)
                fh = logging.FileHandler(logfile)
                ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
                fh.setFormatter(ff)
                self._logger.addHandler(fh)
    
            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)

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
        # first obtain shared lock; may help to avoid race conditions between
        # exclusive locks (REQUIRES THOROUGH TESTING)
        if self._shlock(f):
            fcntl.lockf(f, fcntl.LOCK_EX)
    
        return True

    def _unlock(self, f):
        """Remove exclusive or shared lock on file.

        WARNING: It is very rare that this is necessary, since a file must be unlocked
        before it is closed. Furthermore, locks disappear when a file is closed anyway.
        This method will remain here for now, but may be removed in the future if
        not needed (likely).

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
        """Decorator for opening state file for reading and applying shared lock.
        
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
        for appending and that an exclusive lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

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
    
    @staticmethod
    def _read_data(func):
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
                self.handle = pandas.HDFStore(self.filename, 'r')
                self._shlock(self.handle._handle)
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
            return out

        return inner
    
    @staticmethod
    def _write_data(func):
        """Decorator for opening data file for writing and applying exclusive lock.
        
        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            self.handle = pandas.HDFStore(self.filename, 'a')
            self._exlock(self.handle._handle)
            try:
                out = func(self, *args, **kwargs)
            finally:
                self.handle.close()
            return out


        return inner

class ContainerFile(File):
    """Container file object; syncronized access to Container data.

    """
    class _Meta(tables.IsDescription):
        """Table definition for metadata.

        All strings limited to hardcoded size for now.

        """
        # unique identifier for container
        uuid = tables.StringCol(36)

        # user-given name of container
        name = tables.StringCol(128)

        # container type; Sim or Group
        containertype = tables.StringCol(36)

        # eventually we would like this to be generated dynamically
        # meaning, size of location string is size needed, and meta table
        # is regenerated if any of its strings need to be (or smaller)
        # When Coordinator generates its database, it uses largest string size
        # needed
        location = tables.StringCol(256)

        # version of MDSynthesis object file data corresponds to 
        # allows future-proofing of old objects so that formats of new releases
        # can be automatically built from old ones
        version = tables.StringCol(36)

    class _Coordinator(tables.IsDescription):
        """Table definition for coordinator info.

        This information is kept separate from other metadata to allow the
        Coordinator to simply stack tables to populate its database. It doesn't
        need entries that store its own path.

        Path length fixed size for now.
        """
        # absolute path of coordinator
        abspath = tables.StringCol(256)
        
    class _Tags(tables.IsDescription):
        """Table definition for tags.

        """
        tag = tables.StringCol(36)

    class _Categories(tables.IsDescription):
        """Table definition for categories.

        """
        category = tables.StringCol(36)
        value = tables.StringCol(36)

    def __init__(self, filename, logger=None, **kwargs): 
        """Initialize Container state file.

        This is the base class for all Container state files. It generates 
        data structure elements common to all Containers. It also implements
        low-level I/O functionality.

        :Arguments:
           *filename*
              path to file
           *logger*
              Container's logger instance

        :Keywords:
           *containertype*
              Container type: Sim or Group
           *name*
              user-given name of Container object
           *coordinator*
              directory in which coordinator state file can be found [None]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search

        .. Note:: kwargs passed to :meth:`create`

        """
        super(ContainerFile, self).__init__(filename, logger=logger)

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
           *containertype*
              Container type: Sim or Group
           *name*
              user-given name of Container object
           *coordinator*
              directory in which coordinator state file can be found [None]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search
        """
        containertype = kwargs.pop('containertype', None)

        # metadata table
        self.update_uuid()
        self.update_containertype(containertype)
        self.update_name(kwargs.pop('name', containertype))
        self.update_location()
        self.update_version(kwargs.pop('version', None))

        # coordinator table
        self.update_coordinator(kwargs.pop('coordinator', None))

        # tags table
        tags = kwargs.pop('tags', list())
        self.add_tags(*tags)

        # categories table
        categories = kwargs.pop('categories', dict())
        self.add_categories(**categories)
    
    @File._read_state
    def get_uuid(self):
        """Get Container uuid.
    
        :Returns:
            *uuid*
                unique string for this Container
        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.uuid[0]

    @File._write_state
    def update_uuid(self):
        """Generate new uuid for Container.

        """
        try:
            table = self.handle.get_node('/', 'meta')
            table.cols.uuid[0] = str(uuid4())
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'meta', self._Meta, 'metadata')
            table.row['uuid'] = str(uuid4())
            table.row.append()

    @File._read_state
    def get_name(self):
        """Get Container name.

        :Returns:
            *name*
                name of Container

        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.name[0]

    @File._write_state
    def update_name(self, name):
        """Rename Container.

        :Arugments:
            *name*
                new name of Container

        """
        try:
            table = self.handle.get_node('/', 'meta')
            table.cols.name[0] = name
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'meta', self._Meta, 'metadata')
            table.row['name'] = name
            table.row.append()

    @File._read_state
    def get_containertype(self):
        """Get Container type: Sim or Group.
    
        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.containertype[0]

    @File._write_state
    def update_containertype(self, containertype):
        """Update Container type: Sim or Group.

        Note: will only take 'Sim' or 'Group' as values.

        :Arugments:
            *containertype*
                type of Container: Sim or Group

        """
        if (containertype == 'Sim') or (containertype == 'Group'):
            try:
                table = self.handle.get_node('/', 'meta')
                table.cols.containertype[0] = containertype
            except tables.NoSuchNodeError:
                table = self.handle.create_table('/', 'meta', self._Meta, 'metadata')
                table.row['containertype'] = containertype
                table.row.append()

    @File._read_state
    def get_location(self):
        """Get Container location.

        :Returns:
            *location*
                absolute path to Container directory
    
        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.location[0]

    @File._write_state
    def update_location(self):
        """Update Container location.

        """
        try:
            table = self.handle.get_node('/', 'meta')
            table.cols.location[0] = os.path.dirname(self.filename)
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'meta', self._Meta, 'metadata')
            table.row['location'] = os.path.dirname(self.filename)

    @File._read_state
    def get_version(self):
        """Get Container version.

        :Returns:
            *version*
                version of Container

        """
        table = self.handle.get_node('/', 'meta')
        return table.cols.version[0]

    @File._write_state
    def update_version(self, version):
        """Update version of Container.

        :Arugments:
            *version*
                new version of Container

        """
        try:
            table = self.handle.get_node('/', 'meta')
            table.cols.version[0] = version
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'meta', self._Meta, 'metadata')
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
        return table.cols.abspath[0]

    @File._write_state
    def update_coordinator(self, coordinator):
        """Update Container location.

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
            table = self.handle.create_table('/', 'coordinator', self._Coordinator, 'coordinator information')
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
        return [ x['tag'] for x in table.iterrows() ]
        
    @File._write_state
    def add_tags(self, *tags):
        """Add any number of tags to the Container.

        Tags are individual strings that serve to differentiate Containers from
        one another. Sometimes preferable to categories.

        :Arguments:
           *tags*
              Tags to add. Must be convertable to strings using the str() builtin.

        """
        try:
            table = self.handle.get_node('/', 'tags')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'tags', self._Tags, 'tags')

        # ensure tags are unique (we don't care about order)
        tags = set([ str(tag) for tag in tags ])

        # remove tags already present in metadata from list
        #TODO: more efficient way to do this?
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
        """Delete tags from Container.

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
            tags = set([ str(tag) for tag in tags ])

            # get matching rows
            rowlist = list()
            for row in table:
                for tag in tags:
                    if (row['tag'] == tag):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be removed
            # due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table('/', 'tags', self._Tags, 'tags')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as we
                # delete rows
                for i in rowlist:
                    table.remove_row(i-j)
                    j=j+1

    @File._read_state
    def get_categories(self):
        """Get all categories as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories 
        """
        table = self.handle.get_node('/', 'categories')
        return { x['category']: x['value'] for x in table.iterrows() }

    @File._write_state
    def add_categories(self, **categories):
        """Add any number of categories to the Container.

        Categories are key-value pairs of strings that serve to differentiate
        Containers from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will replace
        the value for that category.

        :Keywords:
            *categories*
                Categories to add. Keyword used as key, value used as value. Both
                must be convertible to strings using the str() builtin.

        """
        try:
            table = self.handle.get_node('/', 'categories')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'categories', self._Categories, 'categories')

        table = self.handle.get_node('/', 'categories')

        # remove categories already present in metadata from dictionary 
        #TODO: more efficient way to do this?
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
        """Delete categories from Container.
    
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
            table = self.handle.create_table('/', 'categories', self._Categories, 'categories')
        else:
            # remove redundant categories from given list if present
            categories = set([ str(category) for category in categories ])

            # get matching rows
            rowlist = list()
            for row in table:
                for category in categories:
                    if (row['category'] == category):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be removed
            # due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table('/', 'categories', self._Categories, 'categories')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as we
                # delete rows
                for i in rowlist:
                    table.remove_row(i-j)
                    j=j+1
    
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

class SimFile(ContainerFile):
    """Main Sim state file.

    This file contains all the information needed to store the state of a
    Sim object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.
    
    """

    class _Topology(tables.IsDescription):
        """Table definition for storing universe topology paths.

        Two versions of the path to a topology are stored: the absolute path
        (abspath) and the relative path from the Sim object's directory
        (relSim). This allows the Sim object to use some heuristically good
        starting points when trying to find missing files using Finder.
        
        """
        abspath = tables.StringCol(255)
        relSim = tables.StringCol(255)

    class _Trajectory(tables.IsDescription):
        """Table definition for storing universe trajectory paths.

        The paths to trajectories used for generating the Universe
        are stored in this table.

        See UniverseTopology for path storage descriptions.

        """
        abspath = tables.StringCol(255)
        relSim = tables.StringCol(255)

    class _Selection(tables.IsDescription):
        """Table definition for storing selections.

        A single table corresponds to a single selection. Each row in the
        column contains a selection string. This allows one to store a list
        of selections so as to preserve selection order, which is often
        required for e.g. structural alignments.

        """
        selection = tables.StringCol(255)

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Sim state file.

        :Arguments:
           *filename*
              path to file
           *logger*
              logger to send warnings and errors to

        """
        super(SimFile, self).__init__(filename, logger=logger, **kwargs)
    
    def create(self, **kwargs):
        """Build Sim data structure.

        :Keywords:
           *name*
              user-given name of Sim object
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
        super(SimFile, self).create(containertype='Sim', **kwargs)

        self._make_universegroup()

    @File._write_state
    def _make_universegroup(self):
        """Make universes group.

        Used only on file creation.

        """
        try:
            group = self.handle.get_node('/', 'universes')
        except tables.NoSuchNodeError:
            group = self.handle.create_group('/', 'universes', 'universes')

    @File._read_state
    def list_universes(self):
        """List universe names.

        :Returns:
            *universes*
                list giving names of all defined universes

        """
        group = self.handle.get_node('/', 'universes')

        return group.__members__

    @File._read_state
    def get_universe(self, universe, path='abspath'):
        """Get topology and trajectory paths for the desired universe.

        :Arguments:
            *universe*
                given name for selecting the universe
            *path*
                type of paths to return; either absolute paths (abspath) or
                paths relative to the Sim object (relSim) ['abspath']

        :Returns:
            *topology*
                path to the topology file
            *trajectory*
                list of paths to trajectory files
                
        """
        # get topology file
        table = self.handle.get_node('/universes/{}'.format(universe), 'topology')
        topology = table.colinstances[path][0]

        # get trajectory files
        table = self.handle.get_node('/universes/{}'.format(universe), 'trajectory')
        trajectory = [ x[path] for x in table.iterrows() ]

        return (topology, trajectory)

    @File._write_state
    def add_universe(self, universe, topology, *trajectory):
        """Add a universe definition to the Sim object.

        A Universe is an MDAnalysis object that gives access to the details
        of a simulation trajectory. A Sim object can contain multiple universe
        definitions (topology and trajectory pairs), since it is often
        convenient to have different post-processed versions of the same
        raw trajectory.

        :Arguments:
            *universe*
                given name for selecting the universe
            *topology*
                path to the topology file
            *trajectory*
                path to the trajectory file; multiple files may be given
                and these will be used in order as frames for the trajectory

        """

        # build this universe's group; if it exists, do nothing 
        try:
            group = self.handle.create_group('/universes', universe, universe, createparents=True)
        except tables.NodeError:
            self.logger.info("Universe definition '{}' already exists. Remove it first.".format(universe))
            return

        # construct topology table 
        table = self.handle.create_table('/universes/{}'.format(universe), 'topology', self._Topology, 'topology')

        # add topology paths to table
        table.row['abspath'] = os.path.abspath(topology)
        table.row['relSim'] = os.path.relpath(topology, self.get_location())
        table.row.append()

        # construct trajectory table
        table = self.handle.create_table('/universes/{}'.format(universe), 'trajectory', self._Trajectory, 'trajectory')

        # add trajectory paths to table
        for segment in trajectory:
            table.row['abspath'] = os.path.abspath(segment)
            table.row['relSim'] = os.path.relpath(segment, self.get_location())
            table.row.append()

        # construct selection group
        group = self.handle.create_group('/universes/{}'.format(universe), 'selections', 'selections')

    @File._write_state
    def del_universe(self, universe):
        """Delete a universe definition.

        Deletes any selections associated with the universe.

        :Arguments:
            *universe*
                name of universe to delete
        """
        self.handle.remove_node('/universes', universe, recursive=True)

    @File._read_state
    def list_selections(self, universe):
        """List selection names.

        :Arguments:
            *universe*
                name of universe the selections apply to

        :Returns:
            *selections*
                list giving names of all defined selections for the given
                universe

        """
        group = self.handle.get_node('/universes/{}'.format(universe), 'selections')

        return group.__members__

    @File._read_state
    def get_selection(self, universe, handle):
        """Get a stored atom selection for the given universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name to use for the selection

        :Returns:
            *selection*
                list of the selection strings making up the atom selection
        """
        table = self.handle.get_node('/universes/{}/selections'.format(universe), handle)
        selection = [ x['selection'] for x in table.iterrows() ]

        return selection

    @File._write_state
    def add_selection(self, universe, handle, *selection):
        """Add an atom selection definition for the named Universe definition.

        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name to use for the selection
            *selection*
                selection string; multiple strings may be given and their
                order will be preserved, which is useful for e.g. structural 
                alignments

        """
        #TODO: add check for existence of selection table
        #TODO: add check for selections as strings; use finally statements
        # to delete table in case of failure
        # construct selection table
        try:
            table = self.handle.create_table('/universes/{}/selections'.format(universe), handle, self._Selection, handle)
        except tables.NoSuchNodeError:
            self.logger.info("Universe definition '{}' does not exist. Add it first.".format(universe))
            return

        # add trajectory paths to table
        for item in selection:
            table.row['selection'] = item
            table.row.append()
    
    @File._write_state
    def del_selection(self, universe, handle):
        """Delete an atom selection from the specified universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name of the selection

        """
        self.handle.remove_node('/universes/{}/selections'.format(universe), handle)

class GroupFile(ContainerFile):
    """Main Group state file.

    This file contains all the information needed to store the state of a
    Group object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.
    
    """

    class _Members(tables.IsDescription):
        """Table definition for the members of the Group.

        Stores for each member its container type, uuid, and two versions of
        the path to the member container: the absolute path (abspath) and the
        relative path from the Group object's directory (relGroup). This allows
        the Group object to use some heuristically good starting points when
        trying to find missing files using Finder.
        
        """
        # unique identifier for container
        uuid = tables.StringCol(36)

        # container type; Sim or Group
        containertype = tables.StringCol(36)

        abspath = tables.StringCol(255)
        relGroup = tables.StringCol(255)

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Group state file.

        :Arguments:
           *filename*
              path to file
           *logger*
              logger to send warnings and errors to

        """
        super(GroupFile, self).__init__(filename, logger=logger, **kwargs)
    
    def create(self, **kwargs):
        """Build Group data structure.

        :Keywords:
           *name*
              user-given name of Sim object
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
        super(GroupFile, self).create(containertype='Group', **kwargs)

        self._make_membertable()

    @File._write_state
    def _make_membertable(self):
        """Make member table.

        Used only on file creation.

        """
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'members', self._Members, 'members')

    @File._write_state
    def add_member(self, uuid, containertype, location):
        """Add a member to the Group.

        If the member is already present, its location will be updated with
        the given location.

        :Arguments:
            *uuid*
                the uuid of the new member
            *containertype*
                the container type of the new member (Sim or Group)
            *location*
                location of the new member in the filesystem
    
        """
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'members', self._Members, 'members')

        # check if uuid already present
        rownum = [ row.nrow for row in table.where("uuid=='{}'".format(uuid)) ]
        if rownum:
            self.logger.info('Member already present. Updating with new location.')
            table.cols.abspath[rownum[0]] = os.path.abspath(location)
            table.cols.relGroup[rownum[0]] = os.path.relpath(location, self.get_location())
        else:
            table.row['uuid'] = uuid
            table.row['containertype'] = containertype
            table.row['abspath'] = os.path.abspath(location)
            table.row['relGroup'] = os.path.relpath(location, self.get_location())
            table.row.append()

    @File._write_state
    def del_member(self, *uuid):
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
            table = self.handle.create_table('/', 'members', self._Members, 'members')
            
        else:
            # remove redundant uuids from given list if present
            uuids = set([ str(tag) for uuid in uuids ])

            # get matching rows
            #TODO: possibly faster to use table.where
            rowlist = list()
            for row in table:
                for uuid in uuids:
                    if (row['uuid'] == uuid):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be removed
            # due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table('/', 'members', self._Members, 'members')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as we
                # delete rows
                for i in rowlist:
                    table.remove_row(i-j)
                    j=j+1

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

        # check if uuid present
        rownum = [ row.nrow for row in table.where("uuid=='{}'".format(uuid)) ]
        if rownum:
            memberinfo = { x: table.colinstances[x][rownum[0]] for x in table.colinstances.keys() }
        else:
            self.logger.info('No such member in this Group.')
            memberinfo = None

        return memberinfo

    @File._read_state
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving uuids of all members, in order

        """
        table = self.handle.get_node('/', 'members')
        return [ x['uuid'] for x in table.iterrows() ]

    @File._read_state
    def get_members_containertype(self):
        """List containertype for each member.

        :Returns:
            *containertypes*
                list giving containertypes of all members, in order

        """
        table = self.handle.get_node('/', 'members')
        return [ x['containertype'] for x in table.iterrows() ]

    @File._read_state
    def get_members_location(self, path='abspath'):
        """List stored location for each member. 

        :Arguments:
            *path*
                type of paths to return; either absolute paths (abspath) or
                paths relative to the Group object (relGroup) ['abspath']

        :Returns:
            *locations*
                list giving locations of all members, in order

        """
        table = self.handle.get_node('/', 'members')
        return [ x[path] for x in table.iterrows() ]

class DatabaseFile(File):
    """Database file object; syncronized access to Database data.

    """

class DataFile(File):
    """Interface to data files.

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
              Container's logger instance

        """
        super(DataFile, self).__init__(filename, logger=logger)

        # open file for the first time to initialize handle
        self.handle = pandas.HDFStore(self.filename, 'r')
        self.handle.close()

    @File._write_data
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
        self.handle.put(key, data, format='table', data_columns=True)

    @File._write_data
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
        self.handle.append(key, data, data_columns=True)

    @File._read_data
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

    @File._write_data
    def del_data(self, key):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete
        
        """
        self.handle.remove(key)
    
    @File._read_data
    def list_data(self):
        """List names of all stored datasets.

        Although the true names start with '\' indicating the root of the
        HDF5 data tree, we wish to abstract this away. We remove the leading
        '\' from the output. This shouldn't cause any problems since the
        leading '\' can be neglected when referring to stored objects by name
        using all of pandas.HDFStore's methods anyway.

        """
        keys = self.handle.keys()
        return [ i.lstrip('/') for i in keys ]
