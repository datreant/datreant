"""
Interface classes for state files.

"""

import os
import warnings
import json

from six import string_types

import datreant
from .core import FileSerial


class MixinJSON(object):
    def _deserialize(self, handle):
        return json.load(handle)

    def _serialize(self, record, handle):
        json.dump(record, handle)


class TreantFile(MixinJSON, FileSerial):
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
        super(FileSerial, self).__init__(filename, logger=logger)

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
                tags to add; must be single numbers, strings, or boolean
                values; tags that are not these types are not added

        """
        # ensure tags are unique (we don't care about order)
        # also they must be of a certain set of types
        tags = set([tag for tag in tags
                    if (isinstance(tag, (int, float, string_types, bool)) or
                        tag is None)])

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
                categories to add; keyword used as key, value used as value;
                values must be single numbers, strings, or boolean values;
                values that are not these types are not added

        """
        for key, value in categories.items():
            if (isinstance(value, (int, float, string_types, bool)) or
                    value is None):
                self._record['categories'][key] = value

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


class GroupFile(TreantFile):
    """Main Group state file.

    This file contains all the information needed to store the state of a
    Group object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.

    """
    # add new paths to include them in member searches
    memberpaths = ['abs', 'rel']
    _fields = ['uuid', 'treanttype', 'abs', 'rel']

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Group state file.

        :Arguments:
           *filename*
              path to file
           *logger*
              logger to send warnings and errors to
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search
        """
        super(GroupFile, self).__init__(filename, logger=logger, **kwargs)

    def _init_record(self):
        super(GroupFile, self)._init_record()
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
                list giving treanttype of each member, in order
        """
        return [member[0] for member in self._record['members']]

    @FileSerial._read
    @FileSerial._pull
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        return [member[1] for member in self._record['members']]

    @FileSerial._read
    @FileSerial._pull
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                list containing all paths to member basedirs, in member order
        """
        return [member[2:] for member in self._record['members']]
