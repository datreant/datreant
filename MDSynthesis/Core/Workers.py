"""
Under-the-hood classes.

"""
import os, sys
import time
import yaml
import cPickle
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

class Utilities(object):
    """Lowest level utilities; contains all methods that are common to every 
       MDSynthesis object.

    """
    def open(self, *args, **kwargs):
        """Open a file for i/o and apply an exclusive lock.
    
        Arguments and keywords are passed directly to the open() builtin.
    
        """
        F = File(*args, **kwargs)
        return F
    
    def makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
class Attributes(object):
    """Class for user-defined attributes.

    """
    def __init__(self):
        pass

