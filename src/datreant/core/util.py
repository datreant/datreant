import os


def makedirs(path, mode=0o777, exist_ok=False):
    """Make directories and all parents necessary.

    :Arguments:
        *path*
            directory path to make
    """
    try:
        os.makedirs(str(path), mode=mode)
    except OSError as e:
        # let slide errors that include directories already existing, but
        # catch others
        if e.errno == 17 and exist_ok:
            pass
        else:
            raise e
