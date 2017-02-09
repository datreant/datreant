import os
from pathlib import Path
import six


def makedirs(path):
    """Make directories and all parents necessary.

    :Arguments:
        *path*
            directory path to make
    """
    try:
        os.makedirs(path)
    except OSError as e:
        # let slide errors that include directories already existing, but
        # catch others
        if e.errno == 17:
            pass
        else:
            raise

            
def touch_me(path):
    Path(path).touch()

    
def relpath(path):
    """Returns *path* on Windows, and relative path elsewhere. """

    if os.name == 'nt':
        return path
    else:
        return os.path.relpath(path)
