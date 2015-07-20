"""
User-level functions for manipulating Treants. Used by included shell
scripts.

"""
from datreant.collections import Bundle


def get(*treants):
    """Get instances of the given Treants.

    :Arguments:
        *treants*
            paths to Treants in the filesystem, or Treant instances,
            to return new Treant instances of

    :Returns:
        *instances*
           Bundle of new instances for each Treant given

    """
    return Bundle(*treants)
