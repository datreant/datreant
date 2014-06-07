"""
The Database object stores information on all Containers it is made aware of.
This centralized storage allows Containers to find each other when necessary;
this is especially important for Groups.

"""
import Core

class Coordinator(Core.Workers.ObjectCore):
    """User-level macrocontroller for all Containers in the database.
    
    """

