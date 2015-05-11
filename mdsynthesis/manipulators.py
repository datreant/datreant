"""
User-level functions for manipulating Containers. Used by included shell
scripts.

"""
import mdsynthesis as mds


def get(*containers):
    """Get instances of the given Containers.

    :Arguments:
        *containers*
            paths to Containers in the filesystem, or Container instances,
            to return new Container instances of

    :Returns:
        *instances*
           Bundle of new instances for each Container given

    """
    return mds.Bundle(*containers)
