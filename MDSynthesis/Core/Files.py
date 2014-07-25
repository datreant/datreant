"""
Classes for datafile syncronization. 

"""

import Aggregators
import Workers
import yaml
import pickle
from uuid import uuid4

class File(object):
    """File object base class. Implements file locking and syncronization.
    
    """
    def __init__(self, filename, reader, writer, datastruct=None, logger=None, **kwargs):
        """Create File instance for interacting with file on disk.

        The File object keeps its own cached version of a data structure
        written to file in sync with the file itself. It does this by
        consistently applying locks to files before reading and writing. 
        
        At all times, reading and modifying the data structure is the same as
        reading and writing to the file. Accessing an element of the data
        structure results in obtaining a shared lock, reading, and releasing
        the lock on the file each time. Modifying an element results in
        obtaining an exclusive lock, reading, writing, and releasing the lock
        on the file. All operations are performed atomically to ensure
        unintended overwrites are avoided.

        :Arguments:
           *filename*
              name of file on disk object synchronizes with
           *reader*
              function used to translate file into data structure;
              must take one argument: file stream to read from
           *writer*
              function used to translate data structure into file;
              must take two arguments: the data structure to write and the 
              file stream to write to
           *datastruct*
              data structure to store; overrides definition in :meth:`create`
              this allows for custom data structures without a pre-defined
              definition 
           *logger*
              logger to send warnings and errors to

        .. Note:: kwargs passed to :meth:`create`

        """
        self.filename = filename
        self.lockname = "{}.lock".format(self.file.name) 

        self.reader = reader
        self.writer = writer

        # if given, use data structure and write to file
        # if none given, check existence of file and read it in if present
        # else create a new data structure and file from definition
        if datastruct:
            self.data = datastruct
            self.locked_write()
        elif self.check_existence():
            self.locked_read()
        else:
            self.create()
            self.locked_write()

    def create(self, **kwargs):
        """Build data structure.

        This is a placeholder function, since each specific File use-case
        will have a different data structure definition.

        """
        self.data = None

    def read(self):
        """Read contents of file into data structure.

        .. Note:: file not locked in this method. Must be done externally.

        :Returns:
           *success*
              True if read successful
        """

        with open(self.filename, 'r') as f:
            self.data = self.reader(f)
        
        return self.compare()

    def write(self):
        """Write data structure to file.
    
        .. Note:: file not locked in this method. Must be done externally.

        :Returns:
           *success*
              True if write successful
        """
        with open(self.filename, 'w') as f:
            self.writer(self.data, f)

        return self.compare()

    def lock(self):
        """Get exclusive lock on file.

        The lock is just a symlink of the target file, since making a symlink
        appears to be an atomic operation on most platforms. This is important,
        since creating a symlink also serves to check if one is already present
        all in one operation.

        :Returns:
           *success*
              True if lockfile successfully created
        """
        # if lockfile already present, wait until it disappears; make lock
        while True:
            try:
                os.symlink(self.filename, self.lockname):
                break
            except OSError:
                time.sleep(1)

        # return lockfile existence
        return os.path.exists(self.lockname)

    def unlock(self):
        """Remove exclusive lock on file.

        Before removing the lock, checks that the data structure is
        the same as what is on file.

        :Returns:
           *success*
              True if comparison successful 
        """
        # check that python representation matches file
        success = self.compare()

        if success:
            os.remove(self.lockname)

        return success

    def lockit(self, func):
        """Decorator for applying a lock around the given method.

        Applying this decorator to a method will ensure that a file lock is
        obtained before that method is executed. It also ensures that the
        lock is removed after the method returns.

        """
        def inner(*args, **kwargs):
            self.lock()
            func(*args, **kwargs)
            self.unlock()

        return inner

    @self.lockit
    def locked_write(self):
        """Lock file, write to file, unlock file.
        
        """
        return self.write()

    @self.lockit
    def locked_read(self):
        """Lock file, write to file, unlock file.

        """
        return self.read()

    def compare(self):
        """Compare data structure with file contents.

        :Returns:
           *same*
              True if file synchronized with data structure
        """
        with open(self.filename, 'r') as f:
            datatemp = self.reader(f)
        
        return self.data == datatemp
    
    def check_existence(self):
        """Check for existence of file.
    
        """
        return os.path.exists(self.filename)

class ContainerFile(File):
    """Container file object; syncronized access to Container data.

    """
    def __init__(self, location, logger, classname, **kwargs): 
        """Initialize Container state file.

        This is the base class for all Container state files. It generates 
        data structure elements common to all Containers.

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

        super(ContainerFile, self).__init__(filename, reader=yaml.load, writer=yaml.dump, logger=logger,
                classname=classname, location=location)

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

class SimFile(ContainerFile):
    """Main Sim state file.

    This file contains all the information needed to store the state of a
    Sim object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well the data structure definition.
    It also defines the format of the file, i.e. the writer and reader
    used to manage it.
    
    """
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

