"""
Classes for datafile syncronization. 

"""

import Aggregators
import Workers

class File(object):
    """File object class. Implements file locking and syncronization.

    """
    def __init__(self, filename):
        """Create File instance for interacting with file on disk.

        The File object keeps its own cached version of a data structure
        written to file in sync with the file itself. It does this by
        consistently applying locks to files before reading and writing. 
        
        At all times, reading and modifying the data structure is the same as
        reading and writing to the file. Accessing an element of the data
        structure results in locking, reading, and unlocking the file each
        time. Modifying an element results in locking, reading, writing, and
        unlocking. All operations are performed atomically to ensure unintended
        overwrites are avoided.

        """
        self.filename = filename
        self.data = self.read(self.filename)

        self.exLockname = "{}.exlock".format(self.file.name) 
        self.shLockname = "{}.shlock".format(self.file.name) 
        self.lock()

    def __enter__(self):
        return self.file

    def __exit__(self, *args):
        self.close()

    def exLock(self):
        """Get exclusive lock on file.

        The lock is just a symlink of the target file, since making a symlink
        appears to be an atomic operation on most platforms. This is important,
        since creating a symlink also serves to check if one is already present
        all in one operation.
    
        """
        # if lockfile already present, wait until it disappears
        while True:
            try:
                os.symlink(self.filename, self.lockname):
                break
            except OSError:
                time.sleep(2)

        # open lockfile (will appear in filesystem)
        self.lockfile = open(self.lockname, 'w')

    def close(self):
        """Close file stream.
        
        """
        self.file.close()
        self.lockfile.close()
        os.remove(self.lockname)

class ContainerFile(File):
    """Container file object; syncronized access to Container data.

    """

class OperatorFile(File):
    """Operator file object; syncronized access to Operator data.

    """

class DatabaseFile(File):
    """Database file object; syncronized access to Database data.

    """

class DataFile(object):
    """Universal datafile interface.

    Allows for safe reading and writing of datafiles, which can be of a wide
    array of formats. Handles the details of conversion from pythonic data
    structure to persistent file form.

    """

