"""Core classes and functions for handling state backends.

"""

import os
import sys
import fcntl
import warnings
import json
from functools import wraps
from contextlib import contextmanager
import six

if six.PY2:
    FileNotFoundError = IOError


class BaseFile(object):
    """Generic File object base class. Implements file locking and reloading
    methods.

    All files in datreant should be accessible by high-level methods
    without having to worry about simultaneous reading and writing by other
    processes. The File object includes methods for ensuring shared and
    exclusive locks are consistently applied before reads and writes,
    respectively. It handles any other low-level tasks for maintaining file
    integrity.

    :Arguments:
        *filename*
            name of file on disk object corresponds to

    Example
    -------

    """

    def __init__(self, filename, **kwargs):
        self.filename = os.path.abspath(filename)
        self.handle = None
        self.fd = None
        self.fdlock = None

        # we apply locks to a proxy file to avoid creating an HDF5 file
        # without an exclusive lock on something; important for multiprocessing
        proxy = "." + os.path.basename(self.filename) + ".proxy"
        self.proxy = os.path.join(os.path.dirname(self.filename), proxy)

        # we create the file if it doesn't exist; if it does, an exception is
        # raised and we catch it; this is necessary to ensure the file exists
        # so we can use it for locks
        try:
            fd = os.open(self.proxy, os.O_CREAT | os.O_EXCL)
            os.close(fd)
        except OSError as e:
            # if we get the error precisely because the file exists, continue
            if e.errno == 17:
                pass
            else:
                raise

    def get_location(self):
        """Get File basedir.

        :Returns:
            *location*
                absolute path to File basedir

        """
        return os.path.dirname(self.filename)

    def delete(self):
        """Delete this file and its proxy file.

        This file instance will be unusable after this operation.

        """
        with self.write():
            os.remove(self.filename)
            os.remove(self.proxy)

    @contextmanager
    def read(self):
        # if we already have any lock, proceed
        if self.fdlock:
            yield self.handle
        else:
            self._apply_shared_lock()

            try:
                # open the file using the actual reader
                self.handle = self._open_file_r()
                yield self.handle
            finally:
                self.handle.close()
                self._release_lock()

    @contextmanager
    def write(self):
        # if we already have an exclusive lock, proceed
        if self.fdlock == 'exclusive':
            yield self.handle
        else:
            self._apply_exclusive_lock()

            # open the file using the actual writer
            self.handle = self._open_file_w()
            try:
                yield self.handle
            finally:
                self.handle.close()
                self._release_lock()

    def exists(self):
        return os.path.exists(self.filename)

    def _apply_shared_lock(self):
        self._open_fd_r()
        self._shlock(self.fd)
        self.fdlock = 'shared'

    def _apply_exclusive_lock(self):
        self._open_fd_rw()
        self._exlock(self.fd)
        self.fdlock = 'exclusive'

    def _release_lock(self):
        self._unlock(self.fd)
        self._close_fd()
        self.fdlock = None

    def _shlock(self, fd):
        """Get shared lock on file.

        Using fcntl.lockf, a shared lock on the file is obtained. If an
        exclusive lock is already held on the file by another process,
        then the method waits until it can obtain the lock.

        :Arguments:
            *fd*
                file descriptor
        """
        fcntl.lockf(fd, fcntl.LOCK_SH)

    def _exlock(self, fd):
        """Get exclusive lock on file.

        Using fcntl.lockf, an exclusive lock on the file is obtained. If a
        shared or exclusive lock is already held on the file by another
        process, then the method waits until it can obtain the lock.

        :Arguments:
            *fd*
                file descriptor
        """
        fcntl.lockf(fd, fcntl.LOCK_EX)

    def _unlock(self, fd):
        """Remove exclusive or shared lock on file.

        WARNING: It is very rare that this is necessary, since a file must be
        unlocked before it is closed. Furthermore, locks disappear when a file
        is closed anyway.  This method will remain here for now, but may be
        removed in the future if not needed (likely).

        :Arguments:
            *fd*
                file descriptor
        """
        fcntl.lockf(fd, fcntl.LOCK_UN)

    def _open_fd_r(self):
        """Open read-only file descriptor for application of advisory locks.

        Because we need an active file descriptor to apply advisory locks to a
        file, and because we need to do this before opening a file with
        the apprpriate interface (e.g. PyTables), we open
        a separate file descriptor to the same file and apply the locks
        to it.

        """
        self.fd = os.open(self.proxy, os.O_RDONLY)

    def _open_fd_rw(self):
        """Open read-write file descriptor for application of advisory locks.

        """
        self.fd = os.open(self.proxy, os.O_RDWR)

    def _close_fd(self):
        """Close file descriptor used for application of advisory locks.

        """
        # close file descriptor for locks
        os.close(self.fd)
        self.fd = None

    def _close(self):
        """Close file.

        Not to be used except for debugging files.

        """
        self.handle.close()
        self._unlock(self.fd)
        self.fdlock = None
        self._close_fd()


class File(BaseFile):
    def _open_file_r(self):
        return open(self.filename)

    def _open_file_w(self):
        return open(self.filename, 'w')


@contextmanager
def atomic_write(fname, file_type=File):
    with file_type(fname).write() as fh:
        yield fh


@contextmanager
def read(fname, file_type=File):
    with file_type(fname).read() as fh:
        yield fh
