"""
Interface classes for state files.

"""

import os
import warnings
import json
from collections import defaultdict

from six import string_types
from six.moves import zip

import datreant
from .core import FileSerial


class MixinJSON(object):
    def _deserialize(self, handle):
        return json.load(handle)

    def _serialize(self, record, handle):
        json.dump(record, handle, indent=4)


class TreantFile(MixinJSON, FileSerial):
    """Treant state file.

    This is the base class for all Treant state files. It generates data
    structure elements common to all Treants. It also implements low-level
    I/O functionality.

    :Arguments:
        *filename*
            path to file
        *logger*
            Treant's logger instance

    :Keywords:
        *treanttype*
            Treant type
        *categories*
            user-given dictionary with custom keys and values; used to
            give distinguishing characteristics to object for search
        *tags*
            user-given list with custom elements; used to give
            distinguishing characteristics to object for search

    .. Note:: kwargs passed to :meth:`create`

    """
    def __init__(self, filename, logger=None, **kwargs):
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
        self._record = {'tags': list(), 'categories': dict()}

    def create(self, **kwargs):
        """Build state file and common data structure elements.

        :Keywords:
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search
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
    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        return self._record['version']

    # TODO: need a proper schema update mechanism
    @FileSerial._write
    def update_schema(self):
        """Update schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            version = self._record['version']
        except KeyError:
            version = datreant.core.__version__

        return version

    @FileSerial._write
    def update_version(self, version):
        """Update version of Treant.

        :Arugments:
            *version*
                new version of Treant
        """
        self._record['version'] = version

    @FileSerial._read
    def get_tags(self):
        """Get all tags as a list.

        :Returns:
            *tags*
                list of all tags
        """
        return self._record['tags']

    @FileSerial._write
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
    def del_tags(self, tags=None, all=False):
        """Delete tags from Treant.

        :Arguments:
            *tags*
                An iterable of tags to delete.
            *all*
                When True, delete all tags [``False``]

        """
        if all:
            self._record['tags'] = list()
        elif tags:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])
            for tag in tags:
                # remove tag; if not present, continue anyway
                try:
                    self._record['tags'].remove(tag)
                except ValueError:
                    pass

    @FileSerial._read
    def get_categories(self):
        """Get all categories as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories
        """
        return self._record['categories']

    @FileSerial._write
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
    def del_categories(self, categories=None, all=False):
        """Delete categories from Treant.

        :Arguments:
            *categories*
                Iterable of category keys to delete.
            *all*
                When True, delete all categories [``False``]

        """
        if all:
            self._record['categories'] = dict()
        elif categories:
            for key in categories:
                # continue even if key not already present
                self._record['categories'].pop(key, None)


class GroupFile(TreantFile):
    """Group state file.

    This is the interface class for Group state files. It generates data
    structure elements for storing information on other Treants as Group
    members.

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

    # add new paths to include them in member searches
    memberpaths = ['abspath', 'relpath']
    _fields = ['uuid', 'treanttype']
    _fields.extend(memberpaths)

    def _init_record(self):
        super(GroupFile, self)._init_record()
        self._record['members'] = list()

    @FileSerial._write
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
        member_rec = {'uuid': uuid,
                      'treanttype': treanttype,
                      'abspath': os.path.abspath(basedir),
                      'relpath': os.path.relpath(
                          basedir, self.get_location())}

        # check if uuid already present
        uuids = [member['uuid'] for member in self._record['members']]

        if uuid in uuids:
            self._record['members'][uuids.index(uuid)] = member_rec

        else:
            # add to record
            self._record['members'].append(member_rec)

    @FileSerial._write
    def del_members(self, uuids=None, all=False):
        """Remove members from the Group.

        :Arguments:
            *uuids*
                An iterable of uuids of the members to remove
            *all*
                When True, remove all members [``False``]

        """
        if all:
            self._record['members'] = list()
        elif uuids:
            # remove redundant uuids from given list if present
            uuids = set([str(uuid) for uuid in uuids])

            # get matching rows
            # TODO: possibly faster to use table.where
            memberlist = list()
            for i, member in enumerate(self._record['members']):
                for uuid in uuids:
                    if (member['uuid'] == uuid):
                        memberlist.append(i)

            memberlist.sort()
            j = 0
            # delete matching entries; have to use j to shift the register as
            # we remove entries
            for i in memberlist:
                self._record['members'].pop(i - j)
                j = j + 1

    @FileSerial._read
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
            if member['uuid'] == uuid:
                memberinfo = member

        return memberinfo

    @FileSerial._read
    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                dict giving full member data, with fields as keys and in member
                order
        """
        out = defaultdict(list)

        for member in self._record['members']:
            for key in self._fields:
                out[key].append(member[key])

        return out

    @FileSerial._read
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving treanttype of each member, in order
        """
        return [member['uuid'] for member in self._record['members']]

    @FileSerial._read
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        return [member['treanttype'] for member in self._record['members']]

    @FileSerial._read
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                list of dicts giving all paths to member basedirs, in member
                order
        """
        return [member.fromkeys(memberpaths)
                for member in self._record['members']]
