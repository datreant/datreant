import os


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


def touch_me(fname, mode=0o666, dir_fd=None, **kwargs):
    if os.name == 'nt':
        flags = os.O_CREAT | os.O_APPEND
        with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as filename:
            os.utime(filename.fileno() if os.utime in os.supports_fd else fname,
                     dir_fd=None if os.supports_fd else dir_fd, **kwargs)
