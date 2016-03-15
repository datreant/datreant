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


def discover(dirpath='.', depth=None, treantdepth=None):
    """Find all Treants within given directory, recursively.

    Parameters
    ----------
    dirpath : string
        Directory within which to search for Treants.
    depth : int
        Maximum directory depth to tolerate while traversing in search of
        Treants. ``None`` indicates no depth limit.
    treantdepth : int
        Maximum depth of Treants to tolerate while traversing in search
        of Treants. ``None`` indicates no Treant depth limit.

    Returns
    -------
    found : Bundle
        Bundle of found Treants.

    """
    from .collections import Bundle
    found = list()

    startdepth = len(dirpath.split(os.sep))
    treantdirs = set()

    for root, dirs, files in scandir.walk(dirpath):
        # depth check; if too deep, next iteration
        if depth and len(root.split(os.sep)) - startdepth > depth:
            continue

        # Treant depth checking
        if treantdepth:

            # remove Treant dirs from our set of them if we've backed out
            for treantdir in list(treantdirs):
                if treantdir not in root:
                    treantdirs.remove(treantdir)

            # actual depth check
            if len(treantdirs) > treantdepth:
                continue

        for treanttype in _TREANTS:
            outnames = fnmatch.filter(files,
                                      "{}.*.json".format(treanttype))

            if treantdepth and outnames:
                treantdirs.add(root)

            paths = [os.path.join(root, file) for file in outnames]
            found.extend(paths)

    return Bundle(found)
