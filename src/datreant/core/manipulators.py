"""
User-level functions for manipulating Treants.

"""
import os
import scandir
import fnmatch

from . import _TREANTS


def get(*treants):
    """Get instances of the given Treants.

    Parameters
    ----------
    treants : str
        Filesystem paths to Treants in the filesystem, or Treant instances,
        to return new Treant instances for.

    Returns
    -------
    instances : Bundle
        Bundle of new instances for each Treant given.

    """
    from .collections import Bundle
    return Bundle(*treants)


def discover(dirpath='.'):
    """Find all Treants within given directory, recursively.

    Parameters
    ----------
    dirpath : string
        Directory within which to search for Treants.

    Returns
    -------
    found : Bundle
        Bundle of found Treants.

    """
    from .collections import Bundle
    found = list()
    for root, dirs, files in scandir.walk(dirpath):
        for treanttype in _TREANTS:
            outnames = fnmatch.filter(files,
                                      "{}.*.json".format(treanttype))
            paths = [os.path.join(root, file) for file in outnames]
            found.extend(paths)

    return Bundle(found)
