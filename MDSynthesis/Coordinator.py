"""
The Coordinator object stores information on all Containers it is made aware of.
This centralized storage allows Containers to find each other more easily when
necessary; this is especially important for Groups. It also provides an interface
for querying the whole set of Containers it is aware of and returning them;
this can remove the need to manually go hunting for Containers in the filesystem.

"""
import Core

class Coordinator(object):
    """Macrocontroller for all Containers in the database.

    """

    def __init__(self, database, **kwargs):
        """Generate a Coordinator from a new or existing database.

        :Arguments:
          *database*
              directory containing a Database file; if no Database file is
              found, a new one is generated

        """

        # initialize aggregators
        #self.add =
        #self.remove = 
        #self.move = 
        pass

    def select(self):
        """
        """
        pass

    
