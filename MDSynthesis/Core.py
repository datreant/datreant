"""
Base classes for :mod:`MDSynthesis` objects.

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

metafile = 'metadata.yaml'
logfile = 'logfile.log'
datafile = 'data.pkl'
dbfile = 'MDSdatabase.yaml'

class ObjectCore(object):
    """Lowest-level mixin; functionality common to all MDSynthesis objects.
    
    """
    def __init__(self):
        """Low-level attribute initialization.

        """
        self.util = Utilities()

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

class Aggregator(object)
    """Core functionality for information aggregators.

    """

class Info(Aggregator):
    """Interface for accessing metadata and status information.

    """
    def __init__(self):
        """
        """
        self.name = 

class Data(Aggregator):
    """Interface for accessing Operator-generated data.

    Combines the results from multiple files, since data for a given Operator
    can be split across many pickled files.

    """

class Bunch(object):
    def __init__(self, odict):
        adict = dict(odict)
        for key in adict:
            if type(adict[key]) is dict:
                adict[key] = RwBunch(adict[key])
        self.__dict__ = adict
