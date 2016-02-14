"""
Limbs are interfaces for accessing stored data, as well as querying
the state of an object.

"""
import os
from six import string_types, with_metaclass
from collections import defaultdict

from . import filesystem
from .collections import Bundle
from . import _TREELIMBS, _LIMBS


class _TreeLimbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _TREELIMBS[limbname] = cls


class _Limbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _LIMBS[limbname] = cls


class TreeLimb(with_metaclass(_TreeLimbmeta, object)):
    """Core functionality for Tree limbs.

    TreeLimbs are meant to attach to Trees, which lack a state file. Since
    Treants are subclasses of Tree, they will inherit attachments of TreeLimbs
    to that class. A TreeLimb, unlike a Limb, does not need access to state
    file data.

    """
    # name used when attached to a Tree's namespace
    _name = 'limb'

    def __init__(self, tree):
        self._tree = tree


class Limb(with_metaclass(_Limbmeta, object)):
    """Core functionality for Treant limbs.

    """
    # name used when attached to a Treant's namespace
    _name = 'limb'

    def __init__(self, treant):
        self._treant = treant

    @property
    def _logger(self):
        return self._treant._logger


class Tags(Limb):
    """Interface to tags.

    """
    _name = 'tags'

    def __init__(self, treant):
        super(Tags, self).__init__(treant)

        # init state if tags not already there;
        # if read-only, check that they are there,
        # and raise exception if they are not
        try:
            with self._treant._write:
                try:
                    self._treant._state['tags']
                except KeyError:
                    self._treant._state['tags'] = list()
        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['tags']
                except KeyError:
                    raise KeyError(
                            ("Missing 'tags' data, and cannot write to "
                             "Treant '{}'".format(self._treant.filepath)))

    def __repr__(self):
        return "<Tags({})>".format(self._list())

    def __str__(self):
        tags = self._list()
        agg = "Tags"
        majsep = "="
        seplength = len(agg)

        if not tags:
            out = "No Tags"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for i in xrange(len(tags)):
                out = out + "'{}'\n".format(tags[i])
        return out

    def __iter__(self):
        return self._list().__iter__()

    def __len__(self):
        return len(self._list())

    def _list(self):
        """Get all tags for the Treant as a list.

        :Returns:
            *tags*
                list of all tags
        """
        with self._treant._read:
            tags = self._treant._state['tags']

        tags.sort()
        return tags

    def add(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
        one another. Sometimes preferable to categories.

        :Arguments:
           *tags*
              Tags to add. Must be convertable to strings using the str()
              builtin.  May also be a list of tags.

        """
        outtags = list()
        for tag in tags:
            if isinstance(tag, list):
                outtags.extend(tag)
            else:
                outtags.append(tag)

        with self._treant._write:
            # ensure tags are unique (we don't care about order)
            # also they must be of a certain set of types
            tags = set([tag for tag in outtags
                        if (isinstance(tag,
                                       (int, float, string_types, bool)) or
                            tag is None)])

            # remove tags already present in metadata from list
            tags = tags.difference(set(self._treant._state['tags']))

            # add new tags

            self._treant._state['tags'].extend(tags)

    def remove(self, *tags):
        """Remove tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.
        """
        with self._treant._write:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])
            for tag in tags:
                # remove tag; if not present, continue anyway
                try:
                    self._treant._state['tags'].remove(tag)
                except ValueError:
                    pass

    def purge(self):
        """Remove all tags from Treant.

        """
        with self._treant._write:
            self._treant._state['tags'] = list()


class Categories(Limb):
    """Interface to categories.

    """
    _name = 'categories'

    def __init__(self, treant):
        super(Categories, self).__init__(treant)

        # init state if categories not already there;
        # if read-only, check that they are there,
        # and raise exception if they are not
        try:
            with self._treant._write:
                try:
                    self._treant._state['categories']
                except KeyError:
                    self._treant._state['categories'] = dict()
        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['categories']
                except KeyError:
                    raise KeyError(
                            ("Missing 'categories' data, and cannot write to "
                             "Treant '{}'".format(self._treant.filepath)))

    def __repr__(self):
        return "<Categories({})>".format(self._dict())

    def __str__(self):
        categories = self._dict()
        agg = "Categories"
        majsep = "="
        seplength = len(agg)

        if not categories:
            out = "No Categories"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for key in categories.keys():
                out = out + "'{}': '{}'\n".format(key, categories[key])
        return out

    def __getitem__(self, key):
        """Get value at given key.

        :Arguments:
            *key*
                key of value to return

        :Returns:
            *value*
                value corresponding to given key
        """
        categories = self._dict()
        return categories[key]

    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        outdict = {key: value}
        self.add(outdict)

    def __delitem__(self, category):
        """Remove category from Treant.

        """
        self.remove(category)

    def __iter__(self):
        return self._dict().__iter__()

    def __len__(self):
        return len(self._dict())

    def _dict(self):
        """Get all categories for the Treant as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories

        """
        with self._treant._read:
            return self._treant._state['categories']

    def add(self, *categorydicts, **categories):
        """Add any number of categories to the Treant.

        Categories are key-value pairs of strings that serve to differentiate
        Treants from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        :Keywords:
            *categorydict*
                dict of categories to add; keys used as keys, values used as
                values. Both keys and values must be convertible to strings
                using the str() builtin.
            *categories*
                Categories to add. Keyword used as key, value used as value.
                Both must be convertible to strings using the str() builtin.

        """
        outcats = dict()
        for categorydict in categorydicts:
            if isinstance(categorydict, dict):
                outcats.update(categorydict)
            else:
                raise TypeError("Invalid arguments; non-keyword" +
                                " arguments must be dicts")

        outcats.update(categories)

        with self._treant._write:
            for key, value in outcats.items():
                if (isinstance(value, (int, float, string_types, bool)) or
                        value is None):
                    self._treant._state['categories'][key] = value

    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        """
        with self._treant._write:
            for key in categories:
                # continue even if key not already present
                self._treant._state['categories'].pop(key, None)

    def purge(self):
        """Remove all categories from Treant.

        """
        with self._treant._write:
            self._treant._state['categories'] = dict()

    def keys(self):
        """Get category keys.

        :Returns:
            *keys*
                keys present among categories
        """
        with self._treant._read:
            return self._treant._state['categories'].keys()

    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among categories
        """
        with self._treant._read:
            return self._treant._state['categories'].values()


class MemberBundle(Limb, Bundle):
    """Member manager for Groups.

    """
    _name = 'members'

    # add new paths to include them in member searches
    _memberpaths = ['abspath', 'relpath']
    _fields = ['uuid', 'treanttype']
    _fields.extend(_memberpaths)

    def __init__(self, treant):
        super(MemberBundle, self).__init__(treant)

        # init state if members not already there;
        # if read-only, check that they are there,
        # and raise exception if they are not
        try:
            with self._treant._write:
                try:
                    self._treant._state['members']
                except KeyError:
                    self._treant._state['members'] = list()
        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['members']
                except KeyError:
                    raise KeyError(
                            ("Missing 'members' data, and cannot write to "
                             "Treant '{}'".format(self._treant.filepath)))

        # member Treant cache
        self._cache = dict()
        self._searchtime = 10

    def __repr__(self):
        return "<MemberBundle({})>".format(self._list())

    def __str__(self):
        names = self.names
        treanttypes = self.treanttypes
        agg = "Members"
        majsep = "="
        seplength = len(agg)

        if not names:
            out = "No Members"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for i, name, treanttype in zip(xrange(len(names)),
                                           names,
                                           treanttypes):
                out = out + "{}\t{}:\t{}\n".format(i, treanttype, name)

        return out

    def _add_member(self, uuid, treanttype, basedir):
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
        with self._treant._write:
            # check if uuid already present
            uuids = [member['uuid'] for member in
                     self._treant._state['members']]

            if uuid not in uuids:
                self._treant._state['members'].append(
                        {'uuid': uuid,
                         'treanttype': treanttype,
                         'abspath': os.path.abspath(basedir),
                         'relpath': os.path.relpath(
                             basedir, self._treant.location)})

    def _del_members(self, uuids=None, all=False):
        """Remove members from the Group.

        :Arguments:
            *uuids*
                An iterable of uuids of the members to remove
            *all*
                When True, remove all members [``False``]

        """
        with self._treant._write:
            if all:
                self._treant._state['members'] = list()
            elif uuids:
                # remove redundant uuids from given list if present
                uuids = set([str(uuid) for uuid in uuids])

                # get matching rows
                # TODO: possibly faster to use table.where
                memberlist = list()
                for i, member in enumerate(self._treant._state['members']):
                    for uuid in uuids:
                        if (member['uuid'] == uuid):
                            memberlist.append(i)

                memberlist.sort()
                j = 0
                # delete matching entries; have to use j to shift the register
                # as we remove entries
                for i in memberlist:
                    self._treant._state['members'].pop(i - j)
                    j = j + 1

    def _get_member(self, uuid):
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
        with self._treant._read:
            for member in self._treant._state['members']:
                if member['uuid'] == uuid:
                    memberinfo = member

        return memberinfo

    def _get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                dict giving full member data, with fields as keys and in member
                order
        """
        out = defaultdict(list)

        with self._treant._read:
            for member in self._treant._state['members']:
                for key in self._fields:
                    out[key].append(member[key])

        return out

    def _get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving treanttype of each member, in order
        """
        with self._treant._read:
            return [member['uuid'] for member in
                    self._treant._state['members']]

    def _get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        with self._treant._read:
            return [member['treanttype'] for member in
                    self._treant._state['members']]

    def _get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                list of dicts giving all paths to member basedirs, in member
                order
        """
        with self._treant._read:
            return [member.fromkeys(_memberpaths)
                    for member in self._treant._state['members']]
