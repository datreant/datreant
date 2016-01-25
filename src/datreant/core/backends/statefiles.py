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

    def _serialize(self, state, handle):
        json.dump(state, handle)


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

    def _init_state(self):
        self._state = {'tags': list(), 'categories': dict()}


    @FileSerial._read
    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        return self._state['version']

    # TODO: need a proper schema update mechanism
    @FileSerial._write
    def update_schema(self):
        """Update schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            version = self._state['version']
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
        self._state['version'] = version


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

    def _init_state(self):
        super(GroupFile, self)._init_state()
        self._state['members'] = list()

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
        # check if uuid already present
        uuids = [member['uuid'] for member in self._state['members']]

        if uuid not in uuids:
            self._state['members'].append(
                    {'uuid': uuid,
                     'treanttype': treanttype,
                     'abspath': os.path.abspath(basedir),
                     'relpath': os.path.relpath(
                         basedir, self.get_location())})

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
            self._state['members'] = list()
        elif uuids:
            # remove redundant uuids from given list if present
            uuids = set([str(uuid) for uuid in uuids])

            # get matching rows
            # TODO: possibly faster to use table.where
            memberlist = list()
            for i, member in enumerate(self._state['members']):
                for uuid in uuids:
                    if (member['uuid'] == uuid):
                        memberlist.append(i)

            memberlist.sort()
            j = 0
            # delete matching entries; have to use j to shift the register as
            # we remove entries
            for i in memberlist:
                self._state['members'].pop(i - j)
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
        for member in self._state['members']:
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

        for member in self._state['members']:
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
        return [member['uuid'] for member in self._state['members']]

    @FileSerial._read
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        return [member['treanttype'] for member in self._state['members']]

    @FileSerial._read
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                list of dicts giving all paths to member basedirs, in member
                order
        """
        return [member.fromkeys(memberpaths)
                for member in self._state['members']]
