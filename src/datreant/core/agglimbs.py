"""
AggLimbs are limbs specifically built for collections, in particular
:class:`datreant.core.collections.Bundle`. They often correspond to Treant
limbs but serve as aggregators over collections of them.

"""
from six import string_types, with_metaclass

from . import filesystem
from . import _AGGTREELIMBS, _AGGLIMBS
from .collections import Bundle


class _AggTreeLimbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _AGGTREELIMBS[limbname] = cls


class _AggLimbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _AGGLIMBS[limbname] = cls


class AggTreeLimb(with_metaclass(_AggTreeLimbmeta, object)):
    """Core functionality for limbs attached to a View.

    """
    _name = 'aggtreelimb'

    def __init__(self, collection):
        self._collection = collection


class AggLimb(with_metaclass(_AggLimbmeta, object)):
    """Core functionality for limbs attached to a collection of Treants.

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
        return "<AggTags({})>".format(self.all)

    def __str__(self):
        tags = self.all
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
        return self.all.__iter__()

    def __len__(self):
        return len(self.all)

    def __getitem__(self, value):
        return [member.tags[value] for member in self._collection]

    @property
    def any(self):
        """List tags present among at least one Treant in collection.

        """
        tags = [set(member.tags) for member in self._collection]
        out = set.union(*tags)

        out = list(out)
        out.sort()

        return out

    @property
    def all(self):
        """List tags present among all Treants in collection.

        """
        tags = [set(member.tags) for member in self._collection]
        out = set.intersection(*tags)

        out = list(out)
        out.sort()

        return out

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

    def clear(self):
        """Remove all tags from each Treant in collection.

        """
        for member in self._collection:
            member.tags.clear()


class AggCategories(AggLimb):
    """Interface to categories.

    """
    _name = 'categories'

    def __init__(self, collection):
        super(AggCategories, self).__init__(collection)

    def __repr__(self):
        return "<AggCategories({})>".format(self._dict())

    # FIX
    def __str__(self):
        categories = self.all()
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

    # FIXX
    def __getitem__(self, keys):
        """Get values for a given key, list of keys, or set of keys.

        If a string is provided for *key*, a list is returned containing the
        values for each member in the collection in member order.

        If a list of keys is provided, a list of lists is returned, where each
        list is the set of values (in member order) corresponding to its
        respective key in the order the keys are given.

        If a set of keys is provided, a dict of lists is returned, where theach
        dict contains the provided keys and the value for each key is a list of
        the corresponding values (in member order) of each member in the
        collection.

        :Arguments:
            *keys*
                key(s) of value(s) to return
        """
        keys_type = type(keys)
        if isinstance(keys_type, str):
            return [member.categories[keys] for member in self._collection]
        elif isinstance(keys_type, list):
            outlist = []
            for key in keys:
                outlist.append([member.categories[keys] for member in self._collection])
            return outlist
        elif isinstance(keys_type, set):
            outdict = {}
            for key in keys:
                outdict[key] = [member.categories[keys] for member in self._collection]
            return outdict
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or dict" +
                            " of strings.")

    ## FIX
    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        outdict = {key: value}
        self.add(outdict)

    ## FIXX
    def __delitem__(self, category):
        """Remove category from each Treant in collection.

        """
        for member in self._collection:
            member.categories.remove(category)

    ## FIX
    def __iter__(self):
        return self._dict().__iter__()

    ## FIX
    def __len__(self):
        return len(self._dict())

    ## FIXX
    def add(self, *categorydicts, **categories):
        """Add any number of categories to each Treant in collection.

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
        for member in self._collection:
            member.categories.add(*categorydicts, **categories)

    ## FIXX
    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        :Arguments:
            *categories*
                Categories to delete.

        """
        for member in self._collection:
            member.categories.remove(*categories)

    ## FIXX
    def clear(self):
        """Remove all categories from Treant.

        """
        for member in self._collection:
            membmer.categories.clear()

    ## FIXX
    def keys(self):
        """Get the keys present among Categories of each Treant in collection.

        :Returns:
            *keys*
                keys present among Categories
        """
        outlist = []
        for member in self._collection:
            outlist.append(member.categories.keys())
        return outlist
        # Should this return a list (with duplicates), a list of lists (of
        # keys for each Treant), or a set?

    ## FIXX
    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among Categories
        """
        outlist = []
        for member in self._collection:
            outlist.append(member.categories.values())
        return outlist
        # Should this return a list (with duplicates), a list of lists (of
        # keys for each Treant), or a set? See keys() method too!

    ## FIX
    def groupby(self, keys):
        """Return groupings of Treants based on Categories.

        :Arguments:
            *keys*
                keys present among Categories

        """

        # given a inputted "keys", return groupings of treants where each group
        # is based on the category values in common among the treants
