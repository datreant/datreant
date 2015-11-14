"""
Interface classes for pytables state files.

"""

import os
import sys
import fcntl
import logging
import warnings
from functools import wraps

import tables
import h5py
import pandas as pd
import numpy as np

import datreant
from datreant.backends.core import File

# number of characters required for uuids
uuidlength = 36

# max length in characters for all paths
pathlength = 511

# max character length of strings used for handles, tags, categories
namelength = 55


class TreantFile(File):
    """Treant file object; syncronized access to Treant data.

    """
    class _Version(tables.IsDescription):
        """Table definition for storing version number of file schema.

        All strings limited to hardcoded size for now.

        """
        # version of datreant file schema corresponds to allows future-proofing
        # of old objects so that formats of new releases can be automatically
        # built from old ones
        version = tables.StringCol(15)

    class _Coordinator(tables.IsDescription):
        """Table definition for coordinator info.

        This information is kept separate from other metadata to allow the
        Coordinator to simply stack tables to populate its database. It doesn't
        need entries that store its own path.

        Path length fixed size for now.
        """
        # absolute path of coordinator
        abspath = tables.StringCol(pathlength)

    class _Tags(tables.IsDescription):
        """Table definition for tags.

        """
        tag = tables.StringCol(namelength)

    class _Categories(tables.IsDescription):
        """Table definition for categories.

        """
        category = tables.StringCol(namelength)
        value = tables.StringCol(namelength)

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
        # filter NaturalNameWarnings from pytables, when they arrive
        warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

        super(TreantFile, self).__init__(filename, logger=logger)

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

    def _open_file_r(self):
        return tables.open_file(self.filename, 'r')

    def _open_file_w(self):
        return tables.open_file(self.filename, 'a')

    @File._write
    def create(self, **kwargs):
        """Build state file and common data structure elements.

        :Keywords:
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

    @File._read
    def get_version(self):
        """Get Treant version.

        :Returns:
            *version*
                version of Treant

        """
        table = self.handle.get_node('/', 'version')
        return table.cols.version[0]

    # TODO: need a proper schema update mechanism
    @File._write
    def update_schema(self):
        """Update schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            table = self.handle.get_node('/', 'version')
            version = table.cols.version[0]
        except tables.NoSuchNodeError:
            version = datreant.__version__

        return version

    @File._write
    def update_version(self, version):
        """Update version of Treant.

        :Arugments:
            *version*
                new version of Treant
        """
        try:
            table = self.handle.get_node('/', 'version')
            table.cols.version[0] = version
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'version', self._Version, 'version')
            table.row['version'] = version
            table.row.append()

    @File._read
    def get_coordinator(self):
        """Get absolute path to Coordinator.

        :Returns:
            *coordinator*
                absolute path to Coordinator directory

        """
        table = self.handle.get_node('/', 'coordinator')
        out = table.cols.abspath[0]

        if out == 'None':
            out = None
        return out

    @File._write
    def update_coordinator(self, coordinator):
        """Update Treant location.

        :Arguments:
            *coordinator*
                absolute path to Coordinator directory
        """
        try:
            table = self.handle.get_node('/', 'coordinator')
            if coordinator:
                table.cols.abspath[0] = os.path.abspath(coordinator)
            else:
                table.cols.abspath[0] = None
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'coordinator', self._Coordinator,
                'coordinator information')
            if coordinator:
                table.row['abspath'] = os.path.abspath(coordinator)
            else:
                table.row['abspath'] = None
            table.row.append()

    @File._read
    def get_tags(self):
        """Get all tags as a list.

        :Returns:
            *tags*
                list of all tags
        """
        table = self.handle.get_node('/', 'tags')
        return [x['tag'] for x in table.read()]

    @File._write
    def add_tags(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
        one another. Sometimes preferable to categories.

        :Arguments:
            *tags*
                Tags to add. Must be convertable to strings using the str()
                builtin.

        """
        try:
            table = self.handle.get_node('/', 'tags')
        except tables.NoSuchNodeError:
            table = self.handle.create_table('/', 'tags', self._Tags, 'tags')

        # ensure tags are unique (we don't care about order)
        tags = set([str(tag) for tag in tags])

        # remove tags already present in metadata from list
        tags = tags.difference(set(table.read()['tag']))

        # add new tags
        for tag in tags:
            table.row['tag'] = tag
            table.row.append()

    @File._write
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
        table = self.handle.get_node('/', 'tags')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table('/', 'tags', self._Tags, 'tags')

        else:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])

            # TODO: improve performance
            # get matching rows
            rowlist = list()
            for row in table:
                for tag in tags:
                    if (row['tag'] == tag):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'tags', self._Tags, 'tags')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1

    @File._read
    def get_categories(self):
        """Get all categories as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories
        """
        table = self.handle.get_node('/', 'categories')
        return {x['category']: x['value'] for x in table.read()}

    @File._write
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
        try:
            table = self.handle.get_node('/', 'categories')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'categories', self._Categories, 'categories')

        table = self.handle.get_node('/', 'categories')

        # remove categories already present in metadata from dictionary
        # TODO: more efficient way to do this?
        for row in table:
            for key in categories.keys():
                if (row['category'] == key):
                    row['value'] = str(categories[key])
                    row.update()
                    # dangerous? or not since we are iterating through
                    # categories.keys() and not categories?
                    categories.pop(key)

        # add new categories
        for key in categories.keys():
            table.row['category'] = key
            table.row['value'] = str(categories[key])
            table.row.append()

    @File._write
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
        table = self.handle.get_node('/', 'categories')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table(
                '/', 'categories', self._Categories, 'categories')
        else:
            # remove redundant categories from given list if present
            categories = set([str(category) for category in categories])

            # get matching rows
            rowlist = list()
            for row in table:
                for category in categories:
                    if (row['category'] == category):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'categories', self._Categories, 'categories')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1


class GroupFile(TreantFile):
    """Main Group state file.

    This file contains all the information needed to store the state of a
    Group object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.

    """
    # add new paths to include them in member searches
    memberpaths = ['abspath', 'relCont']

    class _Members(tables.IsDescription):

        """Table definition for the members of the Group.

        Stores for each member its treant type, uuid, and two versions of
        the path to the member treant: the absolute path (abspath) and the
        relative path from the Group object's directory (relCont). This allows
        the Group object to use some heuristically good starting points when
        trying to find missing files using a Foxhound.

        """
        # unique identifier for treant
        uuid = tables.StringCol(uuidlength)

        # treant type
        treanttype = tables.StringCol(namelength)

        abspath = tables.StringCol(pathlength)
        relCont = tables.StringCol(pathlength)

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
        super(GroupFile, self).__init__(filename, logger=logger, **kwargs)

    def create(self, **kwargs):
        """Build Group data structure.

        :Keywords:
           *coordinator*
              directory in which Coordinator state file can be found [``None``]
           *categories*
              user-given dictionary with custom keys and values; used to
              give distinguishing characteristics to object for search
           *tags*
              user-given list with custom elements; used to give distinguishing
              characteristics to object for search

        .. Note:: kwargs passed to :meth:`create`

        """
        super(GroupFile, self).create(treanttype='Group', **kwargs)

        self._make_membertable()

    @File._write
    def _make_membertable(self):
        """Make member table.

        Used only on file creation.

        """
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

    @File._write
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
        try:
            table = self.handle.get_node('/', 'members')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

        # check if uuid already present
        rownum = [row.nrow for row in table.where("uuid=='{}'".format(uuid))]
        if rownum:
            # if present, update location
            table.cols.abspath[rownum[0]] = os.path.abspath(basedir)
            table.cols.relCont[rownum[0]] = os.path.relpath(
                    basedir, self.get_location())
        else:
            table.row['uuid'] = uuid
            table.row['treanttype'] = treanttype
            table.row['abspath'] = os.path.abspath(basedir)
            table.row['relCont'] = os.path.relpath(
                    basedir, self.get_location())
            table.row.append()

    @File._write
    def del_member(self, *uuid, **kwargs):
        """Remove a member from the Group.

        :Arguments:
            *uuid*
                the uuid(s) of the member(s) to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        table = self.handle.get_node('/', 'members')
        purge = kwargs.pop('all', False)

        if purge:
            table.remove()
            table = self.handle.create_table(
                '/', 'members', self._Members, 'members')

        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uid) for uid in uuid])

            # get matching rows
            # TODO: possibly faster to use table.where
            rowlist = list()
            for row in table:
                for uuid in uuids:
                    if (row['uuid'] == uuid):
                        rowlist.append(row.nrow)

            # must include a separate condition in case all rows will be
            # removed due to a limitation of PyTables
            if len(rowlist) == table.nrows:
                table.remove()
                table = self.handle.create_table(
                    '/', 'members', self._Members, 'members')
            else:
                rowlist.sort()
                j = 0
                # delete matching rows; have to use j to shift the register as
                # we delete rows
                for i in rowlist:
                    table.remove_row(i - j)
                    j = j + 1

    @File._read
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
        table = self.handle.get_node('/', 'members')
        fields = table.dtype.names

        memberinfo = None
        for row in table.where("uuid == '{}'".format(uuid)):
            memberinfo = row.fetch_all_fields()

        if memberinfo:
            memberinfo = {x: y for x, y in zip(fields, memberinfo)}

        return memberinfo

    @File._read
    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                structured array giving full member data, with
                each row corresponding to a member
        """
        table = self.handle.get_node('/', 'members')
        return table.read()

    @File._read
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                array giving treanttype of each member, in order
        """
        table = self.handle.get_node('/', 'members')
        return table.read()['uuid']

    @File._read
    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                array giving treanttype of each member, in order
        """
        table = self.handle.get_node('/', 'members')
        return table.read()['treanttype']

    @File._read
    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                structured array containing all paths to member basedirs
        """
        table = self.handle.get_node('/', 'members')
        return table.read()[self.memberpaths]
