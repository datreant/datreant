"""
Interface classes for state files.

"""

import os
import warnings

from .core import JSONFile


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
