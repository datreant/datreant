"""Core classes and functions for handling state backends.

"""

import os
import sys
import fcntl
import logging
import warnings
from functools import wraps


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

    statefileclass = _TREANTS[treant]._backendclass

    if not statefileclass:
        raise IOError("No known backend type for file '{}'".format(filename))

    return statefileclass(filename, logger=logger, **kwargs)


class File(object):
    """Generic File object base class. Implements file locking and reloading
    methods.

    """

    def __init__(self, filename, logger=None, **kwargs):
        """Create File instance for interacting with file on disk.

        All files in datreant should be accessible by high-level methods
        without having to worry about simultaneous reading and writing by other
        processes. The File object includes methods and infrastructure for
        ensuring shared and exclusive locks are consistently applied before
        reads and writes, respectively. It handles any other low-level tasks
        for maintaining file integrity.

        :Arguments:
            *filename*
                name of file on disk object corresponds to
            *logger*
                logger to send warnings and errors to

        """
        self.filename = os.path.abspath(filename)
        self.handle = None
        self.fd = None
        self.fdlock = None

        self._start_logger(logger)

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

    def _start_logger(self, logger):
        """Start up the logger.

        """
        # delete current logger
        try:
            del self.logger
        except AttributeError:
            pass

        # log to standard out if no logger given
        if not logger:
            self.logger = logging.getLogger(
                '{}'.format(self.__class__.__name__))
            self.logger.setLevel(logging.INFO)

            if not any([isinstance(x, logging.StreamHandler)
                        for x in self.logger.handlers]):
                ch = logging.StreamHandler(sys.stdout)
                cf = logging.Formatter(
                        '%(name)-12s: %(levelname)-8s %(message)s')
                ch.setFormatter(cf)
                self.logger.addHandler(ch)
        else:
            self.logger = logger

    def _shlock(self, fd):
        """Get shared lock on file.

        Using fcntl.lockf, a shared lock on the file is obtained. If an
        exclusive lock is already held on the file by another process,
        then the method waits until it can obtain the lock.

        :Arguments:
            *fd*
                file descriptor

        :Returns:
            *success*
                True if shared lock successfully obtained
        """
        fcntl.lockf(fd, fcntl.LOCK_SH)

        return True

    def _exlock(self, fd):
        """Get exclusive lock on file.

        Using fcntl.lockf, an exclusive lock on the file is obtained. If a
        shared or exclusive lock is already held on the file by another
        process, then the method waits until it can obtain the lock.

        :Arguments:
            *fd*
                file descriptor

        :Returns:
            *success*
                True if exclusive lock successfully obtained
        """
        fcntl.lockf(fd, fcntl.LOCK_EX)

        return True

    def _unlock(self, fd):
        """Remove exclusive or shared lock on file.

        WARNING: It is very rare that this is necessary, since a file must be
        unlocked before it is closed. Furthermore, locks disappear when a file
        is closed anyway.  This method will remain here for now, but may be
        removed in the future if not needed (likely).

        :Arguments:
            *fd*
                file descriptor

        :Returns:
            *success*
                True if lock removed
        """
        fcntl.lockf(fd, fcntl.LOCK_UN)

        return True

    def _open_fd_r(self):
        """Open read-only file descriptor for application of advisory locks.

        Because we need an active file descriptor to apply advisory locks to a
        file, and because we need to do this before opening a file with
        PyTables due to the risk of caching stale state on open, we open
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

    def _apply_shared_lock(self):
        """Apply shared lock.

        """
        self._open_fd_r()
        self._shlock(self.fd)
        self.fdlock = 'shared'

    def _apply_exclusive_lock(self):
        """Apply exclusive lock.

        """
        self._open_fd_rw()
        self._exlock(self.fd)
        self.fdlock = 'exclusive'

    def _release_lock(self):
        """Apply exclusive lock.

        """
        self._unlock(self.fd)
        self._close_fd()
        self.fdlock = None

    @staticmethod
    def _read(func):
        """Decorator for opening state file for reading and applying shared
        lock.

        Applying this decorator to a method will ensure that the file is opened
        for reading and that a shared lock is obtained before that method is
        executed. It also ensures that the lock is removed and the file closed
        after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.fdlock:
                out = func(self, *args, **kwargs)
            else:
                self._apply_shared_lock()

                # open the file using the actual reader
                self.handle = self._open_file_r()
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
                    self._release_lock()
            return out

        return inner

    @staticmethod
    def _write(func):
        """Decorator for opening state file for writing and applying exclusive lock.

        Applying this decorator to a method will ensure that the file is opened
        for appending and that an exclusive lock is obtained before that method
        is executed. It also ensures that the lock is removed and the file
        closed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.fdlock == 'exclusive':
                out = func(self, *args, **kwargs)
            else:
                self._apply_exclusive_lock()

                # open the file using the actual writer
                self.handle = self._open_file_w()
                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self.handle.close()
                    self._release_lock()
            return out

        return inner

    def _open_r(self):
        """Open file with intention to write.

        Not to be used except for debugging files.

        """
        self._open_fd_r()
        self._shlock(self.fd)
        self.fdlock = 'shared'
        self.handle = self._open_file_r()

    def _open_w(self):
        """Open file with intention to write.

        Not to be used except for debugging files.

        """
        self._open_fd_rw()
        self._exlock(self.fd)
        self.fdlock = 'exclusive'
        self.handle = self._open_file_w()

    def _close(self):
        """Close file.

        Not to be used except for debugging files.

        """
        self.handle.close()
        self._unlock(self.fd)
        self.fdlock = None
        self._close_fd()


class FileSerial(File):
    """File object base class for serialization formats, such as JSON.

    """
    def _open_file_r(self):
        return open(self.filename, 'r')

    def _open_file_w(self):
        return open(self.filename, 'w')

    def read_file(self):
        """Return deserialized representation of file.

        """
        self._apply_shared_lock()

        self.handle = self._open_file_r()
        out = self._deserialize(self.handle)
        self.handle.close()

        self._release_lock()

        return out 

    @staticmethod
    def _read(func):
        """Decorator for applying a shared lock on file and reading contents.

        Applying this decorator to a method will ensure a shared lock and the
        latest version of the data is obtained before that method is executed.
        It also ensures that the lock is removed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.fdlock:
                out = func(self, *args, **kwargs)
            else:
                self._apply_shared_lock()
                try:
                    self._pull_state()
                    out = func(self, *args, **kwargs)
                finally:
                    self._release_lock()
            return out

        return inner

    @staticmethod
    def _write(func):
        """Decorator for applying an exclusive lock on file and modifying
        contents.

        Applying this decorator to a method will ensure an exclusive lock and
        the latest version of the data is obtained before that method is
        executed. It also ensures that changes to the data are written to the
        file and the lock removed after the method returns.

        """
        @wraps(func)
        def inner(self, *args, **kwargs):
            if self.fdlock == 'exclusive':
                out = func(self, *args, **kwargs)
            else:
                self._apply_exclusive_lock()
                try:
                    try:
                        self._pull_state()
                    except IOError:
                        self._init_state()
                    out = func(self, *args, **kwargs)
                    self._push_state()
                finally:
                    self._release_lock()
            return out

        return inner

    def _pull_state(self):
        self.handle = self._open_file_r()
        self._state = self._deserialize(self.handle)
        self.handle.close()

    def _deserialize(self, handle):
        """Deserialize full state from open file handle.

        Override with specific code for deserializing stream from *handle*.
        """
        raise NotImplementedError

    def _push_state(self):
        self.handle = self._open_file_w()
        self._serialize(self._state, self.handle)
        self.handle.close()

    def _serialize(self, state, handle):
        """Serialize full state to open file handle.

        Override with specific code for serializing *state* to *handle*.
        """
        raise NotImplementedError
