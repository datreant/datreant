"""
Classes for datafile syncronization. 

"""

import Aggregators
import Workers
import yaml
import pickle
from uuid import uuid4
import tables
import fcntl

class File(object):
    """File object base class. Implements file locking and reloading methods.
    
    """
    def __init__(self, filename, logger=None, **kwargs):
        """Create File instance for interacting with file on disk.

        All files in MDSynthesis should be accessible by high-level
        methods without having to worry about simultaneous reading and writing by
        other processes. The File object includes methods and infrastructure
        for ensuring shared and exclusive locks are consistently applied before
        reads and writes, respectively.

        :Arguments:
           *filename*
              name of file on disk object corresponds to 
           *logger*
              logger to send warnings and errors to

        """
        self.filename = filename
        self.handle = None

    def shlock(self):
        """Get shared lock on file.

        Using fcntl.lockf, a shared lock on the file is obtained. If an
        exclusive lock is already held on the file by another process,
        then the method waits until it can obtain the lock.

        :Returns:
           *success*
              True if shared lock successfully obtained
        """
        fcntl.lockf(self.handle, fcntl.LOCK_SH)

        return True

    def exlock(self):
        """Get exclusive lock on file.
    
        Using fcntl.lockf, an exclusive lock on the file is obtained. If a
        shared or exclusive lock is already held on the file by another
        process, then the method waits until it can obtain the lock.

        :Returns:
           *success*
              True if exclusive lock successfully obtained
        """
        # first obtain shared lock; may help to avoid race conditions between
        # exclusive locks (REQUIRES THOROUGH TESTING)
        if self.shlock():
            fcntl.lockf(self.handle, fcntl.LOCK_EX)
    
        return True

    def unlock(self):
        """Remove exclusive or shared lock on file.

        :Returns:
           *success*
              True if lock removed
        """
        fcntl.lockf(self.handle, fcntl.LOCK_UN)

        return True

    def check_existence(self):
        """Check for existence of file.
    
        """
        return os.path.exists(self.filename)

class ContainerFile(File):
    """Container file object; syncronized access to Container data.

    """
    class Meta(tables.IsDescription):
        """Table definition for metadata.

        All strings limited to hardcoded size for now.

        """
        # unique identifier for container
        uuid = tables.StringCol(36)

        # user-given name of container
        name = tables.StringCol(36)

        # eventually we would like this to be generated dynamically
        # meaning, size of location string is size needed, and meta table
        # is regenerated if any of its strings need to be (or smaller)
        # When Coordinator generates its database, it uses largest string size
        # needed
        location = tables.StringCol(256)

    class Coordinator(tables.IsDescription):
        """Table definition for coordinator info.

        This information is kept separate from other metadata to allow the
        Coordinator to simply stack tables to populate its database. It doesn't
        need entries that store its own path.

        Path length fixed size for now.
        """
        # absolute path of coordinator
        coordinator = tables.StringCol(256)
        
    class Tags(tables.IsDescription):
        """Table definition for tags.

        """
        tag = tables.StringCol(36)

    class Categories(tables.IsDescription):
        """Table definition for categories.

        """
        category = tables.StringCol(36)
        value = tables.StringCol(36)

    def __init__(self, location, logger, classname, **kwargs): 
        """Initialize Container state file.

        This is the base class for all Container state files. It generates 
        data structure elements common to all Containers. It also implements
        low-level I/O functionality.

        :Arguments:
           *location*
              directory that represents the Container
           *logger*
              Container's logger instance
           *classname*
              Container's class name

        :Keywords:
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
           *details*
              user-given string for object notes

        .. Note:: kwargs passed to :meth:`create`

        """

        filename = os.path.join(location, containerfile)

        super(ContainerFile, self).__init__(filename, logger=logger)

    def create(self, **kwargs):
        """Build common data structure elements.

        :Keywords:
           *classname*
              Container's class name
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
           *details*
              user-given string for object notes
        """
        self.data = {}

        self.data['location'] = kwargs.pop('location')
        self.data['coordinator'] = kwargs.pop('coordinator', None)
        self.data['uuid'] = str(uuid4())

        self.data['class'] = kwargs.pop('classname')
        self.data['name'] = kwargs.pop('name', self.data['class'])

        # cumbersome, but if the given categories isn't a dictionary, we fix it
        self.data['categories'] = kwargs.pop('categories', dict())
        if not isinstance(self.data['categories'], dict):
            self.data['categories'] = dict()

        # if given tags isn't a list, we fix it
        self.data['tags'] = kwargs.pop('tags', list())
        if not isinstance(self.data['tags'], list):
            self.data['tags'] = list()

        # if the given details isn't a string, we fix it
        self.data['details'] = kwargs.pop('details', str())
        if not isinstance(self.data['details'], basestring):
            self.data['details'] = str()
    
    def read(self, func):
        """Decorator for opening file for reading and applying shared lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        def inner(*args, **kwargs):
            self.lock()
            func(*args, **kwargs)
            self.unlock()

        return inner

class SimFile(ContainerFile):
    """Main Sim state file.

    This file contains all the information needed to store the state of a
    Sim object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well the data structure definition.
    It also defines the format of the file, i.e. the writer and reader
    used to manage it.
    
    """

    class UniverseTopology(tables.IsDescription):
        """Table definition for storing universe topology paths.

        Three versions of the path to a topology are stored: the absolute path
        (abspath), the relative path from user's home directory (relhome), and the
        relative path from the Sim object's directory (relSim). This allows the
        Sim object to use some heuristically good starting points trying to find
        missing files using Finder.
        
        """
        abspath = tables.StringCol(255)
        relhome = tables.StringCol(255)
        relSim = tables.StringCol(255)

    class UniverseTrajectory(tables.IsDescription):
        """Table definition for storing universe trajectory paths.

        The paths to trajectories used for generating the Universe
        are stored in this table.

        See UniverseTopology for path storage descriptions.

        """
        abspath = tables.StringCol(255)
        relhome = tables.StringCol(255)
        relSim = tables.StringCol(255)

    class Selection(tables.IsDescription):
        """Table definition for storing selections.

        A single table corresponds to a single selection. Each row in the
        column contains a selection string. This allows one to store a list
        of selections so as to preserve selection order, which is often
        required for structural alignments.

        """
        selection = tables.StringCol(255)

    def __init__(self, location, logger, **kwargs):
        """Initialize Sim state file.

        :Arguments:
           *location*
              directory that represents the Container
           *logger*
              logger to send warnings and errors to

        """
        super(SimFile, self).__init__(location, logger=logger, classname='Sim', **kwargs)
    
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
           *details*
              user-given string for object notes

        .. Note:: kwargs passed to :meth:`create`

        """
        super(SimFile, self).create(**kwargs)

        self.data['universes'] = dict()
        self.data['selections'] = dict()

class DatabaseFile(File):
    """Database file object; syncronized access to Database data.

    """

class DataFile(object):
    """Universal datafile interface.

    Allows for safe reading and writing of datafiles, which can be of a wide
    array of formats. Handles the details of conversion from pythonic data
    structure to persistent file form.

    """
    file.flush(fsync=True)

