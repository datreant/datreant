"""
Interface classes for state files.

"""

import os

from .core import JSONFile


def treantfile(filename, **kwargs):
    """Generate or regenerate the appropriate treant file instance from
    filename.

    :Arguments:
        *filename*
            path to state file (existing or to be created), including the
            filename

    :Returns:
        *treantfile*
            treantfile instance attached to the given file

    """
    from .. import _TREANTS

    treanttype = os.path.basename(filename).split(os.extsep)[0]

    try:
        statefileclass = _TREANTS[treanttype]._backendclass
    except KeyError:
        raise IOError("No known treant type for file '{}'".format(filename))

    return statefileclass(filename, **kwargs)


class TreantFile(JSONFile):
    """Treant state file.

    This is the base class for all Treant state files. It generates data
    structure elements common to all Treants. It also implements low-level
    I/O functionality.

    :Arguments:
        *filename*
            path to file

    """
    def _init_state(self):
        self._state = dict()
