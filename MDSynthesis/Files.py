"""
Classes for datafile syncronization. 

"""
from Core import *

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

class MetaFile(File):
    """Metadata file object; syncronized access to metadata files.

    """


class DataFile(File):
    """Datafile file object; syncronized access to datafiles.

    """



