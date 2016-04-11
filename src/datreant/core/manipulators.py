"""
User-level functions for manipulating Treants.

"""
import os
import scandir
import fnmatch
from six.moves import range

from . import _TREANTS


def discover(dirpath='.', depth=None, treantdepth=None):
    """Find all Treants within given directory, recursively.

    Parameters
    ----------
    dirpath : string, Tree
        Directory within which to search for Treants. May also be an existing
        Tree.
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
    from .trees import Tree

    if isinstance(dirpath, Tree):
        if not dirpath.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        dirpath = dirpath.abspath

    found = list()

    startdepth = len(dirpath.split(os.sep))
    treantdirs = set()

    for root, dirs, files in scandir.walk(dirpath):
        for treanttype in _TREANTS:
            outnames = fnmatch.filter(files,
                                      "{}.*.json".format(treanttype))

            if treantdepth is not None and outnames:
                treantdirs.add(root)

            paths = [os.path.join(root, file) for file in outnames]
            found.extend(paths)

        # depth check; if too deep, empty dirs to avoid downward traversal
        if depth is not None and len(root.split(os.sep)) - startdepth >= depth:
            for i in range(len(dirs)):
                dirs.pop()
            continue

        # Treant depth checking
        if treantdepth is not None:

            # remove Treant dirs from our set of them if we've backed out
            for treantdir in list(treantdirs):
                if treantdir not in root:
                    treantdirs.remove(treantdir)

            # actual depth check; if too deep, empty dirs to avoid downward
            # traversal
            if len(treantdirs) > treantdepth:
                for i in range(len(dirs)):
                    dirs.pop()
                continue

    return Bundle(found)
