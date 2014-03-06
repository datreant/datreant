"""
Aggregators are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.).

An Aggregator is designed to be user friendly on its own, but it can be used as
a backend by a Container, too.

"""

class Aggregator(object):
    """Core functionality for information aggregators.

    """

class Info(Aggregator):
    """Interface for accessing metadata and status information.

    """
    def __init__(self):
        """
        """
        self.name = None

class SimInfo(Info):
    """Sim-specific bindings.
    
    """

class GroupInfo(Info):
    """Group-specific bindings.
    
    """

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
