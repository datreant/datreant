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

    **kwargs passed to treant file ``__init__()`` method

    :Returns:
        *treantfile*
            treantfile instance attached to the given file

    """
    from .. import _TREANTS

    treant = None
    basename = os.path.basename(filename)
    for treanttype in _TREANTS:
        if treanttype in basename:
            treant = treanttype
            break

    if not treant:
        raise IOError("No known treant type for file '{}'".format(filename))

    statefileclass = _TREANTS[treant]._backendclass

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
