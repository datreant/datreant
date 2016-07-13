"""
AggLimbs are limbs specifically built for collections, in particular
:class:`datreant.core.collections.Bundle`. They often correspond to Treant
limbs but serve as aggregators over collections of them.

"""
import itertools
import functools
from six import string_types, with_metaclass

from fuzzywuzzy import process

from . import filesystem
from . import _AGGTREELIMBS, _AGGLIMBS
from .collections import Bundle
from .limbs import Tags, Categories


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


@functools.total_ordering
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

    @staticmethod
    def _setter(self, val):
        """Used for constructing the property when attaching this Limb to a class.

        """
        if isinstance(val, (Tags, list, set)):
            val = list(val)
            self.tags.clear()
            self.tags.add(val)
        else:
            raise TypeError("Can only set with tags, a list, or set")

    def __iter__(self):
        return self.all.__iter__()

    def __len__(self):
        return len(self.all)

    def __getitem__(self, value):
        return [member.tags[value] for member in self._collection]

    def __eq__(self, other):
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) == set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __lt__(self, other):
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) < set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __sub__(self, other):
        """Return a set giving the Tags in `a` that are not in `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) - set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __rsub__(self, other):
        """Return a set giving the Tags in `a` that are not in `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(other) - set(self)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __or__(self, other):
        """Return a set giving the union of Tags `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) | set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __ror__(self, other):
        """Return a set giving the union of Tags `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) | set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __and__(self, other):
        """Return a set giving the intersection of Tags `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) & set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __rand__(self, other):
        """Return a set giving the intersection of Tags `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) & set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __xor__(self, other):
        """Return a set giving the symmetric difference of Tags
        `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) ^ set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

    def __rxor__(self, other):
        """Return a set giving the symmetric difference of Tags
        `a` and `b`.

        """
        if isinstance(other, (AggTags, Tags, set)):
            return set(self) ^ set(other)
        else:
            raise TypeError("Operands must be AggTags, Tags, or a set.")

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
        tag : str or list
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

    def filter(self, tag):
        """Filter Treants matching the given tag expression from a Bundle.

        Parameters
        ----------
        tag : str or list
            Tag or tags to filter Treants.

        Returns
        -------
        Bundle
            Bundle of Treants matching the given tag expression.
        """
        matches = self._collection.tags[tag]
        return self._collection[matches]


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

    @staticmethod
    def _setter(self, val):
        """Used for constructing the property when attaching this Limb to a class.

        """
        if isinstance(val, (Categories, dict)):
            val = dict(val)
            self.categories.clear()
            self.categories.add(val)
        else:
            raise TypeError("Can only set with categories or dict")

    def __getitem__(self, keys):
        """Get values for a given key, list of keys, or set of keys.

        If `keys` is a string specifying one key (for a single category),
        return a list of the values among all Treants (in this collection) for
        that category.

        If `keys` is a list of keys, return a list of lists whose order
        corresponds to the order of the elements in `keys`. Each element in
        `keys` is a key specifying a category; each element in the output is
        a list of the values among all Treants (in this collection) for  the
        category specified by the respective key in `keys`.

        If `keys` is a set of keys, return a dict of lists whose keys are the
        same as those provided in `keys`; the value corresponding to each key
        in the output is a list of values among all Treants (in this
        collection) for the category corresponding to that key.

        Parameters
        ----------
        keys : str, list, set
            Valid key(s) of Categories in this collection.

        Returns
        -------
        list, list of lists, dict of lists
            Values for the (single) specified category when `keys` is str.

        """
        if keys is None:
            return None

        members = self._collection
        if isinstance(keys, (int, float, string_types, bool)):
            k = keys
            return [m.categories[k] if k in m.categories else None
                    for m in members]
        elif isinstance(keys, list):
            return [[m.categories[k] if k in m.categories else None
                    for m in members]
                    for k in keys]
        elif isinstance(keys, set):
            return {k: [m.categories[k] if k in m.categories else None
                    for m in members]
                    for k in keys}
        else:
            raise TypeError("Key must be a string, list of strings, or set"
                            " of strings.")

    def __setitem__(self, key, values):
        """Set the value of categories for each Treant in the collection.

        If `values` is not a sequence and is a valid category type (int,
        string_types, bool, float), then it is broadcasted over all members of
        the collection for the category specified by `key`.

        If `values` is a sequence, it must have the same length as the number
        of members in the collection so that, for each member, the value
        assigned to its category (specified by `key`) is the element in
        `values` whose index matches the index of that member in the
        collection.

        Parameters
        ----------
        key : str
            Valid key for the category whose value should be set to *values*.
        values : str, int, float, bool, list, tuple
            Value(s) for the category specified by `key`.
        """
        if values is None:
            return

        members = self._collection
        if isinstance(values, (int, float, string_types, bool)):
            for m in members:
                m.categories.add({key: values})
        elif isinstance(values, (list, tuple)):
            if len(values) != len(members):
                raise ValueError("Values must be a list of the same length as"
                                 " the number of members in the collection.")
            for m, v in zip(members, values):
                m.categories[key] = v

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
        """Get categories present among at least one Treant in collection.

        Returns
        -------
        dict
            All unique Categories among members.
        """
        keys = [set(member.categories.keys()) for member in self._collection]
        keys = set.union(*keys)

        return {k: [m.categories[k] if k in m.categories else None
                for m in self._collection]
                for k in keys}

    @property
    def all(self):
        """Get categories common to all Treants in collection.

        Returns
        -------
        dict
            Categories common to all members.
        """
        keys = [set(member.categories.keys()) for member in self._collection]
        keys = set.intersection(*keys)

        return {k: [m.categories[k] if k in m.categories else None
                for m in self._collection]
                for k in keys}

    def add(self, categorydict=None, **categories):
        """Add any number of categories to each Treant in collection.

        Categories are key-value pairs that serve to differentiate Treants from
        one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        Keys must be strings.

        Values may be ints, floats, strings, or bools. ``None`` as a value
        will not the existing value for the key, if present.

        Parameters
        ----------
        categorydict : dict
            Dict of categories to add; keys used as keys, values used as
            values.
        categories
            Categories to add. Keyword used as key, value used as value.

        """
        for member in self._collection:
            member.categories.add(categorydict, **categories)

    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        Parameters
        ----------
        categories : str
            Categories to delete.
        """
        for member in self._collection:
            member.categories.remove(*categories)

    def clear(self):
        """Remove all categories from all Treants in collection.

        """
        for member in self._collection:
            member.categories.clear()

    def keys(self, scope='all'):
        """Get the keys present among Treants in collection.

        Parameters
        ----------
        scope : {'all', 'any'}
            Keys to return. 'all' will return only keys found within all
            Treants in the collection, while 'any' will return keys found
            within at least one Treant in the collection.

        Returns
        -------
        keys : list
            Present keys.

        """
        keys = [set(member.categories.keys()) for member in self._collection]

        if scope == 'all':
            out = set.intersection(*keys)
        elif scope == 'any':
            out = set.union(*keys)
        else:
            raise ValueError("Scope must be either 'all' or 'any'")

        return list(out)

    def values(self, scope='all'):
        """Get the category values for all Treants in collection.

        Parameters
        ----------
        scope : {'all', 'any'}
            Keys to return. 'all' will return only keys found within all
            Treants in the collection, while 'any' will return keys found
            within at least one Treant in the collection.

        Returns
        -------
        values : list
            A list of values for each Treant in the collection is returned for
            each key within the given `scope`. The value lists are given in
            the same order as the keys from ``AggCategories.keys``.

        """
        keys = self.keys(scope=scope)
        return self[keys]

    def groupby(self, keys):
        """Return groupings of Treants based on values of Categories.

        If a single category is specified by `keys` (`keys` is neither a list
        nor a set of category names), returns a dict of Bundles whose (new)
        keys are the values of the category specified by `keys`; the
        corresponding Bundles are groupings of members in the collection having
        the same category values (for the category specied by `keys`).

        If `keys` is a list of keys, returns a dict of Bundles whose
        (new) keys are tuples of category values. The corresponding Bundles
        contain the members in the collection that have the same set of
        category values (for the categories specified by `keys`); members in
        each Bundle will have all of the category values specified by the tuple
        for that Bundle's key.

        Parameters
        ----------
        keys : str, list
            Valid key(s) of categories in this collection.

        Returns
        -------
        dict
            Bundles of members by category values.
        """
        if keys is None:
            return None

        members = self._collection
        if isinstance(keys, (string_types)):
            catvals = members.categories[keys]
            groupkeys = [v for v in catvals if v is not None]
            groups = {k: Bundle() for k in groupkeys}

            k = keys
            gen = ((m, m.categories[k]) for m in members if
                   k in m.categories and m.categories[k] in groupkeys)
            for m, catval in gen:
                groups[catval].add(m)
        # Note: redundant code in if/elif block can be consolidated in future
        elif isinstance(keys, list):
            catvals = list(zip(*members.categories[keys]))
            groupkeys = [v for v in catvals if None not in v]
            groups = {k: Bundle() for k in groupkeys}

            gen = ((i, m) for i, m in enumerate(members) if
                   set(keys) <= set(m.categories) and
                   tuple([m.categories[k] for k in keys]) in groupkeys)
            for i, m in gen:
                groups[catvals[i]].add(m)

        else:
            raise TypeError("Keys must be a string or a list of"
                            " strings")
        return groups
