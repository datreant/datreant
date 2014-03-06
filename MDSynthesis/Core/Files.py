"""
Classes for datafile syncronization. 

"""

import Aggregators
import Workers

class File(object):
    """File object class. Implements needed file locking and syncronization.

    """
    def __init__(self, *args, **kwargs):
        """Create File instance for loading files and consistently handling locking.

        """
        self.file = open(*args, **kwargs)
        self.lockname = "{}.lock".format(self.file.name) 
        self.lock()

    def __enter__(self):
        return self.file

    def __exit__(self, *args):
        self.close()

    def lock(self):
        """Get lock.
    
        """
        # if lockfile already present, wait until it disappears
        while os.path.exists(self.lockname):
            time.sleep(5)

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

