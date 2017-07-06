"""
Limbs are interfaces for accessing stored data, as well as querying
the state of an object.

"""
import os
import functools
from six import string_types
from collections import defaultdict

from fuzzywuzzy import process

from .collections import Bundle
from .state import JSONFile
from .selectionparser import parse_selection


class Limb(object):
    """Core functionality for Treant limbs.

    """
    def __init__(self, treant):
        self._treant = treant


class MetadataLimb(Limb):

    # TODO; raise meaningful exceptions when file can't be read, written, read
    # and nonexistent, missing, or treantdir missing

    @property
    def _write(self):
        self._statefile = JSONFile(os.path.join(self._treant._treantdir,
                                                self._statefilename),
                                   init_state=self._init_state)

        return self._statefile.write()

    @property
    def _read(self):
        self._statefile = JSONFile(os.path.join(self._treant._treantdir,
                                                self._statefilename),
                                   init_state=self._init_state)

        return self._statefile.read()


@functools.total_ordering
class Tags(MetadataLimb):
    """Interface to tags.

    """
    _name = 'tags'
    _statefilename = 'tags.json'

    @staticmethod
    def _init_state(jsonfile):
        jsonfile._state = []

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

    def __getitem__(self, value):
        # check if we might have a string to parse into a selection object
        if isinstance(value, string_types):
            value = parse_selection(value)
        return self._getselection(value)

    def _getselection(self, value):
        """get tags according to a selection
        """
        with self._read:
            if isinstance(value, list):
                # a list of tags gives only members with ALL the tags
                fits = all([self._getselection(item) for item in value])
            elif isinstance(value, tuple):
                # a tuple of tags gives members with ANY of the tags
                fits = any([self._getselection(item) for item in value])
            if isinstance(value, set):
                # a set of tags gives only members WITHOUT ALL the tags
                # can be used for `not`, basically
                fits = not all([self._getselection(item) for item in value])
            elif isinstance(value, string_types):
                fits = value in self

            return fits

    def __iter__(self):
        return self._list().__iter__()

    def __len__(self):
        return len(self._list())

    def __eq__(self, other):
        if isinstance(other, (Tags, set, list)):
            return set(self) == set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __lt__(self, other):
        if isinstance(other, (Tags, set, list)):
            return set(self) < set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __sub__(self, other):
        """Return a set giving the Tags in `a` that are not in `b`.

        """
        from .agglimbs import AggTags
        if isinstance(other, (AggTags, Tags, set, list)):
            return set(self) - set(other)
        else:
            raise TypeError("Operands must be AggTags, tags, a set, or list.")

    def __rsub__(self, other):
        """Return a set giving the Tags in `a` that are not in `b`.

        """
        from .agglimbs import AggTags  # may not be necessary
        if isinstance(other, (AggTags, Tags, set, list)):
            return set(other) - set(self)
        else:
            raise TypeError("Operands must be AggTags, tags, a set, or list.")

    def __or__(self, other):
        """Return a set giving the union of Tags `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) | set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __ror__(self, other):
        """Return a set giving the union of Tags `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) | set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __and__(self, other):
        """Return a set giving the intersection of Tags `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) & set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __rand__(self, other):
        """Return a set giving the intersection of Tags `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) & set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __xor__(self, other):
        """Return a set giving the symmetric difference of Tags
        `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) ^ set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def __rxor__(self, other):
        """Return a set giving the symmetric difference of Tags
        `a` and `b`.

        """
        if isinstance(other, (Tags, set, list)):
            return set(self) ^ set(other)
        else:
            raise TypeError("Operands must be tags, a set, or list.")

    def _list(self):
        """Get all tags for the Treant as a list.

        :Returns:
            *tags*
                list of all tags
        """
        with self._read:
            tags = self._statefile._state

        tags.sort()
        return tags

    def add(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
        one another. Sometimes preferable to categories.

        Parameters
        ----------
        tags : str or list
            Tags to add. Must be strings or lists of strings.

        """
        outtags = list()
        for tag in tags:
            if isinstance(tag, (list, set, tuple)):
                outtags.extend(tag)
            else:
                outtags.append(tag)

        with self._write:
            # ensure tags are unique (we don't care about order)
            # also they must be strings
            _outtags = []
            for tag in outtags:
                if not isinstance(tag, string_types):
                    raise ValueError("Only string can be added as tags. Tried "
                                     "to add '{}' which is '{}'".format(
                                         tag, type(tag)))
                _outtags.append(tag)
            outtags = set(_outtags)

            # remove tags already present in metadata from list
            outtags = outtags.difference(set(self._statefile._state))

            # add new tags
            self._statefile._state.extend(outtags)

    def remove(self, *tags):
        """Remove tags from Treant.

        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.
        """
        with self._write:
            # remove redundant tags from given list if present
            tags = set([str(tag) for tag in tags])
            for tag in tags:
                # remove tag; if not present, continue anyway
                try:
                    self._statefile._state.remove(tag)
                except ValueError:
                    pass

    def clear(self):
        """Remove all tags from Treant.

        """
        with self._write:
            self._statefile._state = list()

    def fuzzy(self, tag, threshold=80):
        """Get a tuple of existing tags that fuzzily match a given one.

        Parameters
        ----------
        tags : str or list
            Tag or tags to get fuzzy matches for.
        threshold : int
            Lowest match score to return. Setting to 0 will return every tag,
            while setting to 100 will return only exact matches.

        Returns
        -------
        matches : tuple
            Tuple of tags that match.
        """
        if isinstance(tag, string_types):
            tags = [tag]
        else:
            tags = tag

        matches = []

        for tag in tags:
            matches += [i[0] for i in process.extract(tag, self, limit=None)
                        if i[1] > threshold]

        return tuple(matches)


class Categories(MetadataLimb):
    """Interface to categories.

    """
    _name = 'categories'
    _statefilename = 'categories.json'

    @staticmethod
    def _init_state(jsonfile):
        jsonfile._state = {}

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
        """Get values for given `keys`.

        If `keys` is a string, the single value for that string is returned.

        If `keys` is a list of keys, the values for each key are returned in a
        list, in order by the given keys.

        if `keys` is a set of keys, a dict with the keys as keys and values as
        values is returned.

        Parameters
        ----------
        keys : str, list, set
            Key(s) of value to return.

        Returns
        -------
        values : str, int, float, bool, list, or dict
            Value(s) corresponding to given key(s).

        """
        categories = self._dict()

        if isinstance(keys, (int, float, string_types, bool)):
            return categories[keys]
        elif isinstance(keys, list):
            return [categories[key] for key in keys]
        elif isinstance(keys, set):
            return {key: categories[key] for key in keys}
        else:
            raise TypeError("Key must be a string, list of strings, or set"
                            " of strings.")

    def __setitem__(self, key, value):
        """Set value at given key.

        Parameters
        ----------
        key : str
            Key of value to set.
        value : str, int, float, bool
            Value to set for given key.

        """
        outdict = {key: value}
        self.add(outdict)

    def __delitem__(self, category):
        """Remove category from Treant.

        """
        self.remove(category)

    def __eq__(self, other):
        if isinstance(other, (Categories, dict)):
            return dict(self) == dict(other)
        else:
            raise TypeError("Operands must be categories or dicts.")

    def __req__(self, other):
        if isinstance(other, (Categories, dict)):
            return dict(self) == dict(other)
        else:
            raise TypeError("Operands must be categories or dicts.")

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
        with self._read:
            return self._statefile._state

    def add(self, categorydict=None, **categories):
        """Add any number of categories to the Treant.

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
        categories : dict
            Categories to add. Keyword used as key, value used as value.

        """
        outcats = dict()
        if isinstance(categorydict, (dict, Categories)):
            outcats.update(categorydict)
        elif categorydict is None:
            pass
        else:
            raise TypeError("Invalid arguments; non-keyword"
                            " argument must be dict")

        outcats.update(categories)

        with self._write:
            for key, value in outcats.items():
                if not isinstance(key, string_types):
                    raise TypeError("Keys must be strings.")

                if (isinstance(value, (int, float, string_types, bool))):
                    self._statefile._state[key] = value
                elif value is not None:
                    raise TypeError("Values must be ints, floats,"
                                    " strings, or bools.")

    def remove(self, *categories):
        """Remove categories from Treant.

        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.

        Parameters
        ----------
        categories : str
                Categories to delete.

        """
        with self._write:
            for key in categories:
                # continue even if key not already present
                self._statefile._state.pop(key, None)

    def clear(self):
        """Remove all categories from Treant.

        """
        with self._write:
            self._statefile._state = dict()

    def keys(self):
        """Get category keys.

        :Returns:
            *keys*
                keys present among categories
        """
        with self._read:
            return self._statefile._state.keys()

    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among categories
        """
        with self._read:
            return self._statefile._state.values()
