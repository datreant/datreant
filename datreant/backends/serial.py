"""
Base classes for serial state files, such as YAML, JSON.

"""

import os
import sys
import fcntl
import pickle
import logging
import warnings
from functools import wraps

import numpy as np

import datreant
from .core import File

class FileSerial(File):
    def _open_file_r(self):
        return open(self.filename, 'r')

    def _open_file_w(self):
        return open(self.filename, 'w')

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
                self._open_fd_r()
                self._shlock(self.fd)
                self.fdlock = 'shared'

                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self._unlock(self.fd)
                    self._close_fd()
                    self.fdlock = None
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
                self._open_fd_rw()
                self._exlock(self.fd)
                self.fdlock = 'exclusive'

                try:
                    out = func(self, *args, **kwargs)
                finally:
                    self._unlock(self.fd)
                    self.fdlock = None
                    self._close_fd()
            return out

        return inner

    @staticmethod
    def _pull_push(func):
        @wraps(func)
        def inner(self, *args, **kwargs):
            try:
                self._pull_record()
            except IOError:
                self._init_record()
            out = func(self, *args, **kwargs)
            self._push_record()
            return out
        return inner

    @staticmethod
    def _pull(func):
        @wraps(func)
        def inner(self, *args, **kwargs):
            self._pull_record()
            out = func(self, *args, **kwargs)
            return out
        return inner

    def _pull_record(self):
        self.handle = self._open_file_r()
        self._record = self._deserialize(self.handle)
        self.handle.close()

    def _deserialize(self, handle):
        """Deserialize full record from open file handle.
        """
        raise NotImplementedError

    def _push_record(self):
        self.handle = self._open_file_w()
        self._serialize(self._record, self.handle)
        self.handle.close()

    def _serialize(self, record, handle):
        """Serialize full record to open file handle.
        """
        raise NotImplementedError


class TreantFileSerial(FileSerial):
    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Treant state file.

        This is the base class for all Treant state files. It generates data
        structure elements common to all Treants. It also implements
        low-level I/O functionality.

        :Arguments:
            *filename*
                path to file
            *logger*
                Treant's logger instance

        :Keywords:
            *treanttype*
                Treant type
            *name*
                user-given name of Treant object
            *coordinator*
                directory in which coordinator state file can be found [None]
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search
            *version*
                version of datreant file was generated with

        .. Note:: kwargs passed to :meth:`create`

        """
        super(TreantFileSerial, self).__init__(filename, logger=logger)

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

    def _init_record(self):
        self._record = dict()
        self._record['tags'] = list()
        self._record['categories'] = dict()

    def create(self, **kwargs):
        """Build state file and common data structure elements.

        :Keywords:
            *name*
                user-given name of Treant object
            *coordinator*
                directory in which coordinator state file can be found [None]
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search
            *version*
                version of datreant file was generated with
        """
        # update schema and version of file
        version = self.update_schema()
        self.update_version(version)

        # coordinator table
        self.update_coordinator(kwargs.pop('coordinator', None))

        # tags table
        tags = kwargs.pop('tags', list())
        self.add_tags(*tags)

        # categories table
        categories = kwargs.pop('categories', dict())
        self.add_categories(**categories)

    @FileSerial._read
    @FileSerial._pull
    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        return self._record['version']

    # TODO: need a proper schema update mechanism
    @FileSerial._write
    @FileSerial._pull_push
    def update_schema(self):
        """Update schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            version = self._record['version']
        except KeyError:
            version = datreant.__version__

        return version

    @FileSerial._write
    @FileSerial._pull_push
    def update_version(self, version):
        """Update version of Treant.

        :Arugments:
            *version*
                new version of Treant
        """
        self._record['version'] = version

    @FileSerial._read
    @FileSerial._pull
    def get_coordinator(self):
        """Get absolute path to Coordinator.

        :Returns:
            *coordinator*
                absolute path to Coordinator directory

        """
        return self._record['coordinator']

    @FileSerial._write
    @FileSerial._pull_push
    def update_coordinator(self, coordinator):
        """Update Treant location.

        :Arguments:
            *coordinator*
                absolute path to Coordinator directory
        """
        self._record['coordinator'] = coordinator

    @FileSerial._read
    @FileSerial._pull
    def get_tags(self):
        """Get all tags as a list.

        :Returns:
            *tags*
                list of all tags
        """
        return self._record['tags']

    @FileSerial._write
    @FileSerial._pull_push
    def add_tags(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
        one another. Sometimes preferable to categories.

        :Arguments:
            *tags*
                Tags to add. Must be convertable to strings using the str()
                builtin.

        """
        # ensure tags are unique (we don't care about order)
        tags = set([str(tag) for tag in tags])

        # remove tags already present in metadata from list
        tags = tags.difference(set(self._record['tags']))

        # add new tags
        self._record['tags'].extend(tags)

    @FileSerial._write
    @FileSerial._pull_push
    def del_tags(self, *tags, **kwargs):
        """Delete tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.

        :Keywords:
            *all*
                When True, delete all tags [``False``]

        """
        purge = kwargs.pop('all', False)

        if purge:
            self._record['tags'] = list()
        else:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])
            for tag in tags:
                # remove tag; if not present, continue anyway
                try:
                    self._record['tags'].remove(tag)
                except ValueError:
                    pass

    @FileSerial._read
    @FileSerial._pull
    def get_categories(self):
        """Get all categories as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories
        """
        return self._record['categories']

    @FileSerial._write
    @FileSerial._pull_push
    def add_categories(self, **categories):
        """Add any number of categories to the Treant.

        Categories are key-value pairs of strings that serve to differentiate
        Treants from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        :Keywords:
            *categories*
                Categories to add. Keyword used as key, value used as value.
                Both must be convertible to strings using the str() builtin.

        """
        for key in categories.keys():
            self._record['categories'][key] = str(categories[key])

    @FileSerial._write
    @FileSerial._pull_push
    def del_categories(self, *categories, **kwargs):
        """Delete categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        :Keywords:
            *all*
                When True, delete all categories [``False``]

        """
        purge = kwargs.pop('all', False)

        if purge:
            self._record['categories'] = dict()
        else:
            for key in categories:
                # continue even if key not already present
                self._record['categories'].pop(key, None)


class GroupFileSerial(TreantFileSerial):
    """Main Group state file.

    This file contains all the information needed to store the state of a
    Group object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.

    """
    # add new paths to include them in member searches
    memberpaths = ['abspath', 'relCont']
    _fields = ['uuid', 'treanttype', 'abspath', 'relCont']

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Group state file.

        :Arguments:
           *filename*
              path to file
           *logger*
              logger to send warnings and errors to

        :Keywords:
           *coordinator*
              directory in which coordinator state file can be found [None]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search
        """
        super(GroupFileSerial, self).__init__(filename, logger=logger, **kwargs)

    def _init_record(self):
        super(GroupFileSerial, self)._init_record()
        self._record['members'] = list()

    @FileSerial._write
    @FileSerial._pull_push
    def add_member(self, uuid, treanttype, basedir):
        """Add a member to the Group.

        If the member is already present, its basedir paths will be updated
        with the given basedir.

        :Arguments:
            *uuid*
                the uuid of the new member
            *treanttype*
                the treant type of the new member
            *basedir*
                basedir of the new member in the filesystem

        """
        # check if uuid already present
        uuids = [member[0] for member in self._record['members']]

        if uuid not in uuids:
            self._record['members'].append([uuid,
                                            treanttype, 
                                            os.path.abspath(basedir),
                                            os.path.relpath(
                                                basedir, self.get_location())])

    @FileSerial._write
    @FileSerial._pull_push
    def del_member(self, *uuid, **kwargs):
        """Remove a member from the Group.

        :Arguments:
            *uuid*
                the uuid(s) of the member(s) to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        purge = kwargs.pop('all', False)

        if purge:
            self._record['members'] = list()
        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uid) for uid in uuid])

            # get matching rows
            # TODO: possibly faster to use table.where
            memberlist = list()
            for i, member in enumerate(self._record['members']):
                for uuid in uuids:
                    if (member[0] == uuid):
                        memberlist.append(i)

            memberlist.sort()
            j = 0
            # delete matching entries; have to use j to shift the register as
            # we remove entries 
            for i in memberlist:
                self._record['members'].pop(i - j)
                j = j + 1

    @FileSerial._read
    @FileSerial._pull
    def get_member(self, uuid):
        """Get all stored information on the specified member.

        Returns a dictionary whose keys are column names and values the
        corresponding values for the member.

        :Arguments:
            *uuid*
                uuid of the member to retrieve information for

        :Returns:
            *memberinfo*
                a dictionary containing all information stored for the
                specified member
        """
        memberinfo = None
        for member in self._record['members']:
            if member[0] == uuid:
                memberinfo = member

        if memberinfo:
            memberinfo = {x: y for x, y in zip(self._fields, memberinfo)}

        return memberinfo

    @FileSerial._read
    @FileSerial._pull
    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                dict giving full member data, with fields as keys and in member
                order
        """
        out = {key: [] for key in self._fields}

        for member in self._record['members']:
            for i, key in enumerate(self._fields):
                out[key].append(member[i])

        return out

    @FileSerial._read
    @FileSerial._pull
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                array giving treanttype of each member, in order
        """
        return np.array([member[0] for member in self._record['members']])

    @FileSerial._read
    @FileSerial._pull
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                array giving treanttype of each member, in order
        """
        return np.array([member[1] for member in self._record['members']])

    @FileSerial._read
    @FileSerial._pull
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                structured array containing all paths to member basedirs
        """
        return np.array([member[2:] for member in self._record['members']])
