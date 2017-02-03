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

def fullpath(path):
    if not path:
        return None
    if not isinstance(path, six.string_types):
        path = str(path)
    path = os.path.expanduser(path)
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    return path

def isfullpath(path):
    return path == fullpath(path)


def path_leaf(path):
    return os.path.basename(os.path.normpath(path))
