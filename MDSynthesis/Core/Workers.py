"""
Under-the-hood classes. Mostly a grab-bag of needed functionality,
and possibly a bit messy.

"""
import os, sys
import time
import logging
import pdb
import glob
from uuid import uuid4
from multiprocessing import Process

import Aggregators
import Files

class ObjectCore(object):
    """Lowest-level mixin; functionality common to all MDSynthesis user objects.
    
    """
    def __init__(self):
        """Low-level attribute initialization.

        """
        self.util = Utilities()

class Finder(object):
    """A Finder has methods for locating specified objects based on attributes.

    This object is used by Databases to find Containers, and vice-versa.

    """

    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The directory containing
        the first such file found will be returned. None is returned if no such
        files found.

        :Keywords:
            *startdir*
                directory from which to begin upward search; default is
                Container basedir

        :Returns:
            *database*
                directory of located Database; if no Database found, is None
        
        """
        startdir = kwargs.pop('startdir', None)
        
        if not startdir:
            startdir = self.metadata['basedir']

        # search upward for a database
        startdir = os.path.abspath(startdir)
        directory = startdir
        found = False
        
        self._logger.info("Beginning search for database from {}".format(directory))

        while (directory != '/') and (not found):
            directory, tail = os.path.split(directory)
            candidates = glob.glob(os.path.join(directory, self._databasefile))
            
            if candidates:
                self._logger.info("Database candidate located: {}".format(candidates[0]))
                basedir = os.path.dirname(candidates[0])
                db = Database.Database(basedir)
                found = db._handshake()
        
        if not found:
            self._logger.warning("No database found!")
            basedir = None

        return basedir

class Utilities(object):
    """Lowest level utilities; contains all methods that are common to every 
       MDSynthesis object.

    """
    def makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
    def generate_uuid(self):
        """Generate a 'unique' identifier.

        """
        return str(uuid4())

class Attributes(object):
    """Class for user-defined attributes.

    """
    def __init__(self):
        pass

