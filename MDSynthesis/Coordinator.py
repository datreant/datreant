"""
The Database object stores information on all Containers it is made aware of.
This centralized storage allows Containers to find each other when necessary;
this is especially important for Groups.

"""
import Core

class Coordinator(Core.Workers.ObjectCore):
    """User-level macrocontroller for all Containers in the database.

    The Coordinator allows one to 
    
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

    
