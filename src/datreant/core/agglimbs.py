"""
AggLimbs are limbs specifically built for collections, in particular
:class:`datreant.core.collections.Bundle`. They often correspond to Treant
limbs but serve as aggregators over collections of them.

"""
from six import string_types, with_metaclass

from fuzzywuzzy import process

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
        return "<AggTags({})>".format(list(self.all))

    def __str__(self):
        tags = list(self.all)
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
        """Set of tags present among at least one Treant in collection.

        """
        tags = [set(member.tags) for member in self._collection]
        out = set.union(*tags)

        return out

    @property
    def all(self):
        """Set of tags present among all Treants in collection.

        """
        tags = [set(member.tags) for member in self._collection]
        out = set.intersection(*tags)

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

    def fuzzy(self, tag, threshold=80, scope='all'):
        """Get a tuple of existing tags that fuzzily match a given one.

        Parameters
        ----------
        tags : str or list
            Tag or tags to get fuzzy matches for.
        threshold : int
            Lowest match score to return. Setting to 0 will return every tag,
            while setting to 100 will return only exact matches.
        scope : {'all', 'any'}
            Tags to use. 'all' will use only tags found within all Treants in
            collection, while 'any' will use tags found within at least one
            Treant in collection.

        Returns
        -------
        matches : tuple
            Tuple of tags that match.
        """
        if isinstance(tag, string_types):
            tags = [tag]
        else:
            tags = tag

        if scope == 'all':
            choices = self.all
        elif scope == 'any':
            choices = self.any
        else:
            raise ValueError("Scope can only be 'any' or 'all'")

        matches = []

        for tag in tags:
            matches += [i[0] for i in process.extract(tag, choices, limit=None)
                        if i[1] > threshold]

        return tuple(matches)


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
        collection = self._collection
        if isinstance(keys, str):
            k = keys
            return [m.categories[k] if k in m.categories else None
                    for m in collection]
        elif isinstance(keys, list):
            return [[m.categories[k] if k in m.categories else None
                    for m in collection]
                    for k in keys]
        elif isinstance(keys, set):
            return {k: [m.categories[k] if k in m.categories else None
                    for m in collection]
                    for k in keys}
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
        list of list
            Keys for each member are returned as a list, where the list of the
            lists for each member has the same order as the members in the
            collection.
        """
        return [member.categories.keys() for member in self._collection]

    def values(self):
        """Get the unique category values of all Treants in collection.

        Returns
        -------
        list of list
            Values for each member are returned as a list, where the list of
            the lists for each member has the same order as the members in the
            collection.
        """
        return [member.categories.values() for member in self._collection]

    def groupby(self, keys):
        """Return groupings of Treants based on categories.

        Groups Treants according to a given key, list of keys, or set of keys.
        Mirrors the behavior of ``__getitem()__``, where the output type is
        controlled by the input type of `keys`:

        If `keys` is a string specifying one key (for a single category),
        returns a list of Treants that have that category.

        If `keys` is a list of keys, returns a list of lists whose order
        corresponds to the order of the elements in `keys`. Each element in
        `keys` is a key specifying a category; each element in the output is a
        list of the Treants (in this collection) that have the category
        specified by that key.

        If `keys` is a set of keys, returns a dict of lists whose keys are the
        same as those provided in *keys*; the value corresponding to each key
        in the output is a list of Treants (in this collection) that have the
        Category corresponding to that key.

        Parameters
        ----------
        keys : str, list of str, set of str
            Key(s) of Categories in this collection.

        Returns
        -------
        list of Treants, list of list of Treants, dict of list of Treant
            Treants with the specified (single) category when `keys` is str.

            Groupings of Treants, each grouping a list of Treants, where the
            first grouping contains Treants with the Category specified by the
            first value in `keys`, the second grouping contains Treants for the
            second value in `keys`, etc. when `keys` is a list of str.

            Values in the dict corresponding to each of the provided `keys`
            is a grouping (list) of Treants that have the Category specified by
            that key when `keys` is a set of str.
        """
        collection = self._collection
        if isinstance(keys, str):
            k = keys
            return [m for m in collection if k in m.categories]
        elif isinstance(keys, list):
            return [[m for m in collection if k in m.categories]
                    for k in keys]
        elif isinstance(keys, set):
            return {k: [m for m in collection if k in m.categories]
                    for k in keys}
        else:
            raise TypeError("Invalid argument; argument must be" +
                            " a string, list of strings, or set" +
                            " of strings.")
