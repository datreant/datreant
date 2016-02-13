"""
Interface classes for state files.

"""

import os
import warnings
import json
from collections import defaultdict

import datreant
from .core import JSONFile, FileSerial


def treantfile(filename, logger=None, **kwargs):
    """Generate or regenerate the appropriate treant file instance from
    filename.

    :Arguments:
        *filename*
            path to state file (existing or to be created), including the
            filename
        *logger*
            logger instance to pass to treant file instance

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

    return TreantFile(filename, **kwargs)


class TreantFile(JSONFile):
    """Treant state file.

    This is the base class for all Treant state files. It generates data
    structure elements common to all Treants. It also implements low-level
    I/O functionality.

    :Arguments:
        *filename*
            path to file
        *logger*
            Treant's logger instance

    :Keywords:
        *categories*
            user-given dictionary with custom keys and values; used to
            give distinguishing characteristics to object for search
        *tags*
            user-given list with custom elements; used to give
            distinguishing characteristics to object for search

    .. Note:: kwargs passed to :meth:`create`

    """

    def __init__(self, filename, logger=None, **kwargs):
        super(TreantFile, self).__init__(filename, logger=logger)

        # if file does not exist, it is created; if it does exist, it is
        # updated
        try:
            self.create(**kwargs)
        except OSError:
            # in case the file exists but is read-only; we can't update but may
            # still want to use it
            if os.path.exists(self.filename):
                pass
            # if the file doesn't exist, we still want an exception
            else:
                raise

    def create(self, **kwargs):
        """Build state file and common data structure elements.

        """
        # update schema and version of file
        version = self.update_schema()
        self.update_version(version)

    def _init_state(self):
        self._state = dict()

    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        with self.read():
            return self._state['version']

    # TODO: need a proper schema update mechanism
    def update_schema(self):
        """Update schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        with self.write():
            try:
                version = self._state['version']
            except KeyError:
                version = datreant.core.__version__

            return version

    def update_version(self, version):
        """Update version of Treant.

        :Arugments:
            *version*
                new version of Treant
        """
        with self.write():
            self._state['version'] = version
