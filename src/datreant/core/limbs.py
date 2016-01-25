"""
Limbs are interfaces for accessing stored data, as well as querying
the state of an object.

"""
import six

from . import filesystem
from . import collections
from . import _LIMBS


class _Limbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _LIMBS[limbname] = cls


class Limb(six.with_metaclass(_Limbmeta, object)):
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
            with self._treant._write():
                try:
                    self._treant._state['tags']
                except KeyError:
                    self._treant._state['tags'] = list()
        except (IOError, OSError):
            with self._treant.read():
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
        with self._read():
            tags = self._state['tags']

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

        with self._write():
            # ensure tags are unique (we don't care about order)
            # also they must be of a certain set of types
            tags = set([tag for tag in outtags
                        if (isinstance(tag, (int, float, string_types, bool)) or
                            tag is None)])

            # remove tags already present in metadata from list
            tags = tags.difference(set(self._record['tags']))

            # add new tags

            self._state['tags'].extend(tags)

    def remove(self, *tags):
        """Remove tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.
        """
        with self._write():
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])
            for tag in tags:
                # remove tag; if not present, continue anyway
                try:
                    self._record['tags'].remove(tag)
                except ValueError:
                    pass

    def purge(self):
        """Remove all tags from Treant.

        """
        with self._write():
            self._state['tags'] = list()


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
            with self._treant._write():
                try:
                    self._treant._state['categories']
                except KeyError:
                    self._treant._state['categories'] = dict()
        except (IOError, OSError):
            with self._treant.read():
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
        self._backend.del_categories((category,))

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
        with self._read():
            return self._state['categories']

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

        with self._write():
            for key, value in outcats.items():
                if (isinstance(value, (int, float, string_types, bool)) or
                        value is None):
                    self._state['categories'][key] = value

    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        """
        self._backend.del_categories(categories)

        with self._write():
            for key in categories:
                # continue even if key not already present
                self._state['categories'].pop(key, None)

    def purge(self):
        """Remove all categories from Treant.

        """
        with self._write():
            self._state['categories'] = dict()

    def keys(self):
        """Get category keys.

        :Returns:
            *keys*
                keys present among categories
        """
        with self._read():
            return self._state['categories'].keys()

    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among categories
        """
        with self._read():
            return self._state['categories'].values()


class Members(Limb, collections.CollectionBase):
    """Member manager for Groups.

    """
    _name = 'members'

    def __init__(self, treant):
        super(Members, self).__init__(treant)
        self._cache = dict()

    def __repr__(self):
        return "<Members({})>".format(self._list())

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
