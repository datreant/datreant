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
            for key, value in categories.items():
                out = out + "'{}': '{}'\n".format(key, value)
        return out

    def __getitem__(self, keys):
        """Get values for a given key, list of keys, or set of keys.

        If *keys* is a string specifying one key (for a single category),
        return a list of the values among all Treants (in this collection) that
        have that category.

        If *keys* is a list of keys, return a list of lists whose order
        corresponds to the order of the elements in *keys*. Each element in
        *keys* is a key specifying a category; each element in the output is
        a list of the values among all Treants (in this collection) that have
        the category specified by the respective key in *keys*.

        If *keys* is a set of keys, return a dict of lists whose keys are the
        same as those provided in *keys*; the value corresponding to each key
        in the output is a list of values among all Treants (in this
        collection) that have the category corresponding to that key.

        Parameters
        ----------
        keys : str, list of str, set of str
            Key(s) of Categories in this collection.

        Returns
        -------
        list, list of list, dict of list
            Values for the (single) specified category when *keys* is str.

            Groupings of values, each grouping a list, where the first grouping
            contains the values for all members of the collection corresponding
            to the first value in *keys*, the second grouping contains values
            for all members of the collection corresponding the second value in
            *keys*, etc. when *keys* is a list of str.

            Values in the dict corresponding to each of the provided *keys*
            is a grouping (list) of Treants that have the Category specified by
            that key when *keys* is a set of str.
        """
        keys_type = type(keys)
        if isinstance(keys_type, str):
            return [memb.categories[key] for memb in self._collection]
        elif isinstance(keys_type, list):
            return [[memb.categories[key] for memb in self._collection]
                    for key in keys]
        elif isinstance(keys_type, set):
            return {key: [memb.categories[key] for memb in self._collection]
                    for key in keys}
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or set" +
                            " of strings.")

    def __setitem__(self, key, value):
        """Set the value of a category for all Treants in the collection.

        Parameters
        ----------
        key : str
            Key for the category whose value should be set to *value*.
        value
            The value of the category specified by *key*
        """
        for member in self._collection:
            member.categories.add({key: value})

    def __delitem__(self, category):
        """Remove *category* from each Treant in collection.

        """
        for member in self._collection:
            member.categories.remove(category)

    def __iter__(self):
        """Iterator over Categories common to all Treants in collection.

        """
        return self.all.__iter__()

    def __len__(self):
        """The number of Categories common to all Treants in collection.

        """
        return len(self.all)

    @property
    def any(self):
        """Get Categories present among at least one Treant in collection.

        Returns
        -------
        dict
            All unique Categories among members.
        """
        keys = [set(member.categories.keys()) for member in self._collection]
        everykey = set.union(*keys)

        return {k: [m.categories[k] if k in m.categories else None
                for m in self._collection]
                for k in everykey}

    @property
    def all(self):
        """Get Categories common to all Treants in collection.

        Returns
        -------
        dict
            Categories common to all members.
        """
        keys = [set(member.categories.keys()) for member in self._collection]
        sharedkeys = set.intersection(*keys)

        return {k: [m.categories[k] if k in m.categories else None
                for m in self._collection]
                for k in sharedkeys}

    def add(self, *categorydicts, **categories):
        """Add any number of categories to each Treant in collection.

        Categories are key-value pairs of strings that serve to differentiate
        Treants from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        Parameters
        ----------
        *categorydict
            dict of categories to add; keys used as keys, values used as
            values. Both keys and values must be convertible to strings
            using the str() builtin.
        **categories
            Categories to add. Keyword used as key, value used as value.
            Both must be convertible to strings using the str() builtin.
        """
        for member in self._collection:
            member.categories.add(*categorydicts, **categories)

    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        Parameters
        ----------
        *categories
            Categories to delete.
        """
        for member in self._collection:
            member.categories.remove(*categories)

    def clear(self):
        """Remove all categories from Treant.

        """
        for member in self._collection:
            member.categories.clear()

    def keys(self):
        """Get the unique Categories (keys) of all Treants in collection.

        Returns
        -------
        set
            All unique keys present among all Categories in collection.
        """
        return {member.categories.keys() for member in self._collection}

    def values(self):
        """Get the unique category values of all Treants in collection.

        Returns
        -------
        set
            All unique values present among all Categories in collection.
        """
        return {member.categories.values() for member in self._collection}

    def groupby(self, keys):
        """Return groupings of Treants based on Categories.

        Groups Treants according to a given key, list of keys, or set of keys.
        Mirrors the behavior of *__getitem()__*, where the output type is
        controlled by the input type of *keys*:

        If *keys* is a string specifying one key (for a single category),
        *groupby()* returns a list of Treants that have that category.

        If *keys* is a list of keys, *groupby()* returns a list of lists whose
        order corresponds to the order of the elements in *keys*. Each element
        in *keys* is a key specifying a category; each element in the output is
        a list of the Treants (in this collection) that have the category
        specified by that key.

        If *keys* is a set of keys, *groupby()* returns a dict of lists whose
        keys are the same as those provided in *keys*; the value corresponding
        to each key in the output is a list of Treants (in this collection)
        that have the Category corresponding to that key.

        Parameters
        ----------
        keys : str, list of str, set of str
            Key(s) of Categories in this collection.

        Returns
        -------
        list of Treant, list of list of Treant, dict of list of Treant
            Treants with the specified (single) category when *keys* is str.

            Groupings of Treants, each grouping a list of Treants, where the
            first grouping contains Treants with the Category specified by the
            first value in *keys*, the second grouping contains Treants for the
            second value in *keys*, etc. when *keys* is a list of str.

            Values in the dict corresponding to each of the provided *keys*
            is a grouping (list) of Treants that have the Category specified by
            that key when *keys* is a set of str.
        """
        keys_type = type(keys)
        if isinstance(keys_type, str):
            return [m for m in self._collection
                    if m.categories[key] is not None]
        elif isinstance(keys_type, list):
            return [[m for m in self._collection
                    if m.categories[key] is not None]
                    for key in keys]
        elif isinstance(keys_type, set):
            return {key: [m for m in self._collection
                    if m.categories[key] is not None]
                    for key in keys}
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or set" +
                            " of strings.")
