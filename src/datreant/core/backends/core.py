"""Core classes and functions for handling state backends.

"""

import os
import sys
import fcntl
import warnings
import json
from functools import wraps
from contextlib import contextmanager


class File(object):
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

    def delete(self):
        """Delete this file and its proxy file.

        This file instance will be unusable after this operation.

        """
        with self.write():
            os.remove(self.filename)
            os.remove(self.proxy)


class FileSerial(File):
    """File object base class for serialization formats, such as JSON.

    """
    @property
    def _writebuffer(self):
        wbuffer = ".{}.buffer".format(os.path.basename(self.filename))
        return os.path.join(os.path.dirname(self.filename), wbuffer)

    def _open_file_r(self):
        return open(self.filename, 'r')

    def _open_file_w(self):
        return open(self._writebuffer, 'w')

    def read_file(self):
        """Return deserialized representation of file.

        """
        self._apply_shared_lock()

        self.handle = self._open_file_r()
        out = self._deserialize(self.handle)
        self.handle.close()

        self._release_lock()

        return out

    @contextmanager
    def read(self):
        # if we already have any lock, proceed
        if self.fdlock:
            yield self._state
        else:
            self._apply_shared_lock()
            try:
                self._pull_state()
                yield self._state
            finally:
                self._release_lock()

    @contextmanager
    def write(self):
        # if we already have an exclusive lock, proceed
        if self.fdlock == 'exclusive':
            yield self._state
        else:
            self._apply_exclusive_lock()
            try:
                self._pull_state()
            except IOError:
                self._init_state()
            try:
                yield self._state
                self._push_state()
            finally:
                self._release_lock()

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
        os.rename(self._writebuffer, self.filename)

    def _serialize(self, state, handle):
        """Serialize full state to open file handle.

        Override with specific code for serializing *state* to *handle*.
        """
        raise NotImplementedError


class JSONFile(FileSerial):
    def _deserialize(self, handle):
        return json.load(handle)

    def _serialize(self, state, handle):
        json.dump(state, handle)
