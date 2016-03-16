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
        return "<AggCategories({})>".format(self.all)

    def __str__(self):
        categories = self.all
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

        If a set of keys is provided, a dict of lists is returned, where each
        dict contains the provided keys and the value for each key is a list of
        the corresponding values (in member order) of each member in the
        collection.

        :Arguments:
            *keys*
                key(s) of value(s) to return
        """
        keys_type = type(keys)
        if isinstance(keys_type, str):
            return [memb.categories[key] for memb in self._collection]
        elif isinstance(keys_type, list):
            return [[memb.categories[key] for memb in self._collection]
                    for key in keys]
        elif isinstance(keys_type, set):
            # return {key: [memb.categories[key] for memb in self._collection]
            #         for key in keys}
            # Python 2.6 and earlier
            return dict(
                        (key,
                            [m.categories[key] for m in self._collection])
                        for key in keys)
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or dict" +
                            " of strings.")

    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        for member in self._collection:
            member.categories.add({key: value})

    def __delitem__(self, category):
        """Remove category from each Treant in collection.

        """
        for member in self._collection:
            member.categories.remove(category)

    def __iter__(self):
        """Return an iterator across each unique Category in collection.

        """
        return self.all.__iter__()

    def __len__(self):
        """Return number of unique Categories in collection.

        """
        return len(self.all)

    @property
    def any(self):
        """List Categories present among at least one Treant in collection.

        """
        cats = [set(member.categories) for member in self._collection]
        out = set.union(*cats)

        out = list(out)
        out.sort()

        return out

    @property
    def all(self):
        """List Categories present among all Treants in collection.

        """
        cats = [set(member.categories) for member in self._collection]
        out = set.intersection(*cats)

        out = list(out)
        out.sort()

        return out

    # FIXX
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

    def clear(self):
        """Remove all categories from Treant.

        """
        for member in self._collection:
            member.categories.clear()

    # FIXX - this now has analogous behavior as any() for the keys
    def keys(self):
        """Get the unique Categories (keys) of all Treants in collection.

        :Returns:
            *keys*
                set of keys present among Categories
        """
        out = set()
        for member in self._collection:
            out.add(member.categories.keys())
        return out
        # return {member.categories.keys() for member in self._collection}

    # FIXX - this now has analogous behavior as any() for the values
    def values(self):
        """Get the unique category values of all Treants in collection.

        :Returns:
            *values*
                set of values present among Categories
        """
        out = set()
        for member in self._collection:
            out.add(member.categories.values())
        return out
        # return {member.categories.values() for member in self._collection}

    # FIXX
    def groupby(self, keys):
        """Return groupings of Treants based on Categories.

        Groups Treants according to a given key, list of keys, or set of keys.
        Mirrors the behavior of *__getitem()__*:
            * a string as a key returns a list of Treants with that Category
            * a list of keys returns a list of lists of Treants, where the Nth
              list of Treants are those having the Category defined by the Nth
              key
            * a set of keys returns a dict of list of Treants (analogous to
              inputting a list of keys) where for each key is a list of Treants
              having the Category corresponding to that key

        :Arguments:
            *keys*
                keys present among Categories

        """
        keys_type = type(keys)
        if isinstance(keys_type, str):
            group = []
            for member in self._collection:
                if member.categories[key] is not None:
                    group.append(member)
            return group
        elif isinstance(keys_type, list):
            out = []
            for key in keys:
                group = []
                for member in self._collection:
                    if member.categories[key] is not None:
                        group.append(member)
                out.append(group)
            return out
            # return [[member.categories[key] for member in self._collection]
            #         for key in keys]
        elif isinstance(keys_type, set):
            out = {}
            for key in keys:
                group = []
                for member in self._collection:
                    if member.categories[key] is not None:
                        group.append(member)
                out[key] = group
            return out
            # return {
            #      key: [member.categories[key] for member in self._collection]
            #      for key in keys}
            # Python 2.6 and earlier
            # return dict(
            #     (key,[member.categories[key] for member in self._collection])
            #     for key in keys)
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or dict" +
                            " of strings.")
