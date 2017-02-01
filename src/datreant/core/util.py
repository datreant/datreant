import os
import time

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
    if os.name == 'nt':
        now = time.time()
        try:
            # assume it's there
            os.utime(path, (now, now))
        except os.error:
            # if it isn't, try creating the directory,
            # a file with that name
            os.makedirs(os.path.dirname(path))
            open(path, "w").close()
            os.utime(path, (now, now))
