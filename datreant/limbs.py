"""
Limbs are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Treant) in
ways that are user-friendly.

In short, an Limb is designed to be user friendly on its own, but it can
be used as a backend by a Treant, too.

"""

from . import filesystem
from . import collections


class Limb(object):
    """Core functionality for information limbs.

    """
    # name used when attached to a Treant's namespace
    _name = 'limb'

    def __init__(self, treant):
        self._treant = treant
        self._placeholders()

    def _placeholders(self):
        """Initialize any hidden elements.

        """
        pass

    @property
    def _backend(self):
        return self._treant._backend

    @property
    def _logger(self):
        return self._treant._logger


class Tags(Limb):
    """Interface to tags.

    """
    _name = 'tags'

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
        return self._backend.get_tags().__iter__()

    def __len__(self):
        return len(self._backend.get_tags())

    def _list(self):
        """Get all tags for the Treant as a list.

        :Returns:
            *tags*
                list of all tags
        """
        tags = self._backend.get_tags()
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
        self._backend.add_tags(*outtags)

    def remove(self, *tags, **kwargs):
        """Remove tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.

        :Keywords:
            *all*
                When True, delete all tags [``False``]
        """
        self._backend.del_tags(*tags, **kwargs)


class Categories(Limb):
    """Interface to categories.

    """
    _name = 'categories'

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
        categories = self._backend.get_categories()
        return categories[key]

    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        outdict = {key: value}
        self._backend.add_categories(**outdict)

    def __delitem__(self, category):
        """Remove category from Treant.

        """
        self._backend.del_categories(category)

    def __iter__(self):
        return self._backend.get_categories().__iter__()

    def __len__(self):
        return len(self._backend.get_categories())

    def _dict(self):
        """Get all categories for the Treant as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories

        """
        return self._backend.get_categories()

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
        self._backend.add_categories(**outcats)

    def remove(self, *categories, **kwargs):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        :Keywords:
            *all*
                When True, delete all categories [``False``]

        """
        self._backend.del_categories(*categories, **kwargs)

    def keys(self):
        """Get category keys.

        :Returns:
            *keys*
                keys present among categories
        """
        return self._backend.get_categories().keys()

    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among categories
        """
        return self._backend.get_categories().values()


class Members(Limb, collections.CollectionBase):
    """Member manager for Groups.

    """
    _name = 'members'

    def _placeholders(self):
        # member cache
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


class MemberAgg(object):
    """Core functionality for limbs attached to the Members limb.

    """

    def __init__(self, members):
        self._members = members
