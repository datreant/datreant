"""
AggLimbs are limbs specifically built for collections, in particular
:class:`datreant.core.collections.Bundle`. They often correspond to Treant
limbs but serve as aggregators over collections of them.

"""
import six

from . import filesystem
from . import collections
from . import _AGGLIMBS


class _AggLimbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _AGGLIMBS[limbname] = cls


class AggLimb(six.with_metaclass(_AggLimbmeta, object)):
    """Core functionality for limbs attached to a collection.

    """
    _name = 'agglimb'

    def __init__(self, collection):
        self._collection = collection


class AggTags(AggLimb):
    """Interface to aggregated tags.

    """
    _name = 'tags'

    def __init__(self, collection):
        super(AggTags, self).__init__(collection)

    def __repr__(self):
        return "<AggTags({})>".format(self.list())

    def __str__(self):
        tags = self.list()
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
        return self.list().__iter__()

    def __len__(self):
        return len(self.list())

    def _list(self):
        """Get all tags for the Treants in the collection as a list.

        :Returns:
            *tags*
                list of all tags
        """
        with self._treant._read:
            tags = self._treant._state['tags']

        tags.sort()
        return tags

    def list(self, mode='any'):
        """List tags present among Treants in collection.

        :Arguments:
            *mode*
                'any' returns a list of tags present in at least one member
                of the collection; 'all' returns only tags present in all
                members.

        :Returns:
            *tags*
                list of tags matching `mode`

        """
        tags = [set(member.tags) for member in self._collection]
        if mode == 'any':
            out = set.union(*tags)
        elif mode == 'all':
            out = set.intersection(*tags)

        out = list(out)
        out.sort()

        return out

    @property
    def any(self):
        """List tags present among at least one Treant in collection.

        """
        return self.list('any')

    @property
    def all(self):
        """List tags present among all Treants in collection.

        """
        return self.list('all')

    def add(self, *tags):
        """Add any number of tags to each Treant in collection.

        :Arguments:
           *tags*
              Tags to add. Must be strings or lists of strings.

        """
        for member in self._collection:
            member.tags.add(*tags)

    def remove(self, *tags):
        """Remove tags from each Treant in collection.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.
        """
        for member in self._collection:
            member.tags.remove(*tags)

    def purge(self):
        """Remove all tags from each Treant in collection.

        """
        for member in self._collection:
            member.tags.purge()


class AggCategories(AggLimb):
    """Interface to categories.

    """
    _name = 'categories'

    def __init__(self, collection):
        super(AggCategories, self).__init__(collection)

    def __repr__(self):
        return "<AggCategories({})>".format(self._dict())

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
