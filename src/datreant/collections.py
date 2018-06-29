"""
The Bundle object is the primary manipulator for Treants in aggregate. They are
returned as queries to Bundles. They offer convenience methods
for dealing with many Treants at once. Views give the same kind of aggregation
conveniences for Trees and Leaves.

"""

from __future__ import absolute_import

import os
import functools
from collections import defaultdict

import multiprocessing as mp
import glob
import fnmatch

from six import string_types
from six.moves import zip

from .trees import Tree, Leaf
from .names import TREANTDIR_NAME
from .exceptions import NotATreantError
from .metadata import Tags, Categories, AggTags, AggCategories


@functools.total_ordering
class CollectionMixin(object):
    """Mixin class for collections.

    """

    def __len__(self):
        return len(self._list())

    def __iter__(self):
        return self._list().__iter__()

    def __eq__(self, other):
        try:
            return set(self) == set(other)
        except AttributeError:
            return NotImplemented

    def __lt__(self, other):
        try:
            return set(self) < set(other)
        except AttributeError:
            return NotImplemented

    def __getitem__(self, index):
        memberlist = self._list()

        # if (is a list of bools) OR (a numpy array of bool dtype)
        if ((isinstance(index, list) and
             all(isinstance(item, bool) for item in index)) or
                (hasattr(index, 'dtype') and index.dtype == 'bool')):
            # boolean indexing, either with list or np array
            out = self.__class__([memberlist[i]
                                  for i, val in enumerate(index) if val])
        # if is list or array of ints
        elif (isinstance(index, list) or
              (hasattr(index, 'dtype') and index.dtype == 'int')):
            # fancy indexing, either with list or np array
            out = self.__class__([memberlist[item] for item in index])
        elif isinstance(index, int):
            # an index gets the member at that position
            out = memberlist[index]
        else:
            raise IndexError("Cannot index {} with given values"
                             "".format(self.__class__.__name__))
        return out

    def _membertrees(self):
        return View([member for member in self if isinstance(member, Tree)])

    def leaves(self, hidden=False):
        """Return a View of the files within the member Trees.

        Parameters
        ----------
        hidden : bool
            If True, include hidden files.

        Returns
        -------
        View
            A View giving the files in the member Trees.

        """
        return View([member.leaves(hidden=hidden)
                     for member in self._membertrees()])

    def trees(self, hidden=False):
        """Return a View of the directories within the member Trees.

        Parameters
        ----------
        hidden : bool
            If True, include hidden directories.

        Returns
        -------
        View
            A View giving the directories in the member Trees.

        """
        return View([member.trees(hidden=hidden)
                     for member in self._membertrees()])

    def children(self, hidden=False):
        """Return a View of all files and directories within the member Trees.

        Parameters
        ----------
        hidden : bool
            If True, include hidden files and directories.

        Returns
        -------
        View
            A View giving the files and directories in the member Trees.

        """
        return View([member.children(hidden=hidden)
                     for member in self._membertrees()])

    def glob(self, pattern):
        """Return a View of all child Leaves and Trees of members matching
        given globbing pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match files and directories with

        """
        return View([member.glob(pattern) for member in self
                     if isinstance(member, Tree)])

    def draw(self, depth=None, hidden=False):
        """Print an ASCII-fied visual of all member Trees.

        Parameters
        ----------
        depth : int
            Maximum directory depth to display. ``None`` indicates no limit.
        hidden : bool
            If False, do not show hidden files; hidden directories are still
            shown if they contain non-hidden files or directories.

        """
        for member in self:
            if isinstance(member, Tree):
                member.draw(depth=depth, hidden=hidden)

    def parents(self):
        """Return a View of the parent directories for each member.

        Because a View functions as an ordered set, and some members of this
        collection may share a parent, the View of parents may contain fewer
        elements than this collection.

        """
        return View([member.parent for member in self])

    @property
    def loc(self):
        """Get a View giving Tree/Leaf at `path` relative to each Tree in
        collection.

        Use with getitem syntax, e.g. ``.loc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If directory/file does not exist at the given path, then whether a Tree
        or Leaf is given is determined by the path semantics, i.e. a trailing
        separator ("/").

        """
        if not hasattr(self, "_loc"):
            self._loc = _Loc(self)

        return self._loc

    @property
    def treeloc(self):
        """Get a View giving Tree at `path` relative to each Tree in
        collection.

        Use with getitem syntax, e.g. ``.loc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If the given path resolves to an existing file for any Tree, then a
        ``ValueError`` will be raised.

        """
        if not hasattr(self, "_treeloc"):
            self._treeloc = _TreeLoc(self)

        return self._treeloc

    @property
    def leafloc(self):
        """Get a View giving Leaf at `path` relative to each Tree in
        collection.

        Use with getitem syntax, e.g. ``.loc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If the given path resolves to an existing directory for any Tree, then
        a ``ValueError`` will be raised.

        """
        if not hasattr(self, "_leafloc"):
            self._leafloc = _LeafLoc(self)

        return self._leafloc


class View(CollectionMixin):
    """An ordered set of Trees and Leaves.

    Parameters
    ----------
    vegs : Tree, Leaf, or list
        Trees and/or Leaves to be added, which may be nested lists of Trees
        and Leaves. Trees and Leaves can be given as either objects or
        paths.

    """

    def __init__(self, *vegs):
        self._state = list()
        self._add(*vegs)

    def __repr__(self):
        names = [i.name if isinstance(i, Leaf) else i.name + '/'
                 for i in self._list()]
        return "<View({})>".format(names)

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        A single integer will yield a single member. Lists of either will yield
        a View with members in order by the given items. Giving a basename will
        always yield a View, since more than one Tree or Leaf may have the same
        basename and so they are not guaranteed to be unique.

        A boolean index by way of a list or numpy array can also be used to
        select out members.

        """
        memberlist = self._list()
        # we can take lists of indices or names; these return a
        # View; repeats already not respected since View functions as a
        # set
        if isinstance(index, string_types):
            # a name can be used for indexing
            # always returns a View
            out = View([memberlist[i] for i, name
                        in enumerate(self.names) if name == index])
        elif isinstance(index, slice):
            # we also take slices, obviously
            out = View(*memberlist[index])
        else:
            out = super(View, self).__getitem__(index)
        return out

    def __str__(self):
        out = "<- View ->\n"

        for member in self._list():
            out += "  {}\n".format(member.relpath)

        out += "<- ---- ->"
        return out

    def __add__(self, other):
        """Addition of View with a View or Tree/Leaf View.

        """
        if isinstance(other, (Tree, Leaf, View)):
            return View(self, other)
        else:
            raise TypeError("Right operand must be a Tree, Leaf, or View.")

    def __sub__(self, other):
        """Return a View giving the Treants in `self` that are not in `other`.

        Subtracting a Tree/Leaf from a View also works.

        """
        if isinstance(other, View):
            return View(list(set(self) - set(other)))
        elif isinstance(other, (Tree, Leaf)):
            return View(list(set(self) - set([other])))
        else:
            raise TypeError("Right operand must be a Tree, Leaf, or View.")

    def __or__(self, other):
        """Return a View giving the union of Views `self` and `other`.

        """
        if isinstance(other, View):
            return View(self, other)
        else:
            raise TypeError("Operands must be Views.")

    def __and__(self, other):
        """Return a View giving the intersection of Views `self` and `other`.

        """
        if isinstance(other, View):
            return View(list(set(self) & set(other)))
        else:
            raise TypeError("Operands must be Views.")

    def __xor__(self, other):
        """Return a View giving the symmetric difference of Views `self` and
        `other`.

        """
        if isinstance(other, View):
            return View(list(set(self) ^ set(other)))
        else:
            raise TypeError("Operands must be Views.")

    def _add(self, *vegs):
        """Add any number of members to this collection.

        :Arguments:
            *vegs*
                Trees or Leaves to add; lists, tuples, or other
                Views with Trees or Leaves will also work; strings
                giving a path (existing or not) also work, since these
                are what define Trees and Leaves

        """
        from .trees import Veg, Leaf, Tree
        from .treants import Treant

        outconts = list()
        for veg in vegs:
            if veg is None:
                pass
            elif isinstance(veg, (list, tuple)):
                self._add(*veg)
            elif isinstance(veg, (View, Bundle)):
                self._add(*list(veg))
            elif isinstance(veg, Treant):
                outconts.append(Tree(veg))
            elif isinstance(veg, Veg):
                outconts.append(veg)
            elif (isinstance(veg, string_types) and
                    (os.path.isdir(veg) or veg.endswith(os.sep))):
                tre = Tree(veg)
                outconts.append(tre)
            elif isinstance(veg, string_types):
                tre = Leaf(veg)
                outconts.append(tre)
            else:
                raise TypeError("'{}' not a valid input "
                                "for View".format(veg))

        self._add_members(*outconts)

    def _add_members(self, *members):
        """Add many members at once.

        :Arguments:
            *members*
                list of Trees and Leaves

        """
        for member in members:
            self._add_member(member)

    def _add_member(self, member):
        """Add a member to the View.

        :Arguments:
            *member*
                Tree or Leaf to add

        """
        if member not in self._state:
            self._state.append(member)

    def _list(self):
        """Return a list of members.

        """
        return list(self._state)

    @property
    def names(self):
        """List the basenames for the members in this View.

        """
        return [member.name for member in self]

    @property
    def membertrees(self):
        """A View giving only members that are Trees (or subclasses).

        """
        return self._membertrees()

    @property
    def memberleaves(self):
        """A View giving only members that are Leaves (or subclasses).

        """
        return View([member for member in self if isinstance(member, Leaf)])

    @property
    def abspaths(self):
        """List of absolute paths for the members in this View.

        """
        return [member.abspath for member in self]

    @property
    def relpaths(self):
        """List of relative paths from the current working directory for the
        members in this View.

        """
        return [member.relpath for member in self]

    @property
    def exists(self):
        """List giving existence of each member as a boolean.

        """
        return [member.exists for member in self]

    def map(self, function, processes=1, **kwargs):
        """Apply a function to each member, perhaps in parallel.

        A pool of processes is created for `processes` > 1; for example,
        with 40 members and ``processes=4``, 4 processes will be created,
        each working on a single member at any given time. When each process
        completes work on a member, it grabs another, until no members remain.

        `kwargs` are passed to the given function when applied to each member

        Parameters
        ----------
        function : function
            Function to apply to each member. Must take only a single Treant
            instance as input, but may take any number of keyword arguments.
        processes : int
            How many processes to use. If 1, applies function to each member in
            member order in serial.

        Returns
        -------
        results : list
            List giving the result of the function for each member, in member
            order. If the function returns ``None`` for each member, then only
            ``None`` is returned instead of a list.
        """
        if processes > 1:
            pool = mp.Pool(processes=processes)
            results = dict()
            results = {member.abspath: pool.apply_async(
                    function, args=(member,), kwds=kwargs) for member in self}

            pool.close()
            pool.join()

            output = {key: results[key].get() for key in results}

            # sort by member order
            results = [output[abspath] for abspath in self.abspaths]
        else:
            results = [function(member, **kwargs) for member in self]

        # check if list is all ``None``: if so, we return ``None``
        if all([(i is None) for i in results]):
            results = None

        return results

    def globfilter(self, pattern):
        """Return a View of members that match by name the given globbing
        pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match member names with

        """
        return View([self[name] for name in
                     fnmatch.filter(self.names, pattern)])

    def make(self):
        """Make the Trees and Leaves in this View if they don't already exist.

        Returns
        -------
        View
            This View.

        """
        for member in self:
            member.make()

        return self


class Bundle(CollectionMixin):
    """An ordered set of Treants.

    Parameters
    ----------
    treants : Treant, list
        Treants to be added, which may be nested lists of Treants. Treants
        can be given as either objects or paths to directories that contain
        Treant statefiles. Glob patterns are also allowed, and all found
        Treants will be added to the collection.

    """

    def __init__(self, *treants):
        self._cache = dict()
        self._state = list()

        # add metadata objects
        self._tags = AggTags(self)
        self._categories = AggCategories(self)

        self._add(*treants)

    def __repr__(self):
        return "<Bundle({})>".format(self.names)

    def __str__(self):
        out = "<- Bundle ->\n"

        for member in self._list():
            out += "  {}\n".format(member)

        out += "<- ---- ->"
        return out

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        A single integer will yield a single Treant. Lists of
        either will yield a Bundle with members in order by the given items.
        Giving a name will always yield a Bundle, since names are not
        guaranteed to be unique.

        A boolean index by way of a list or numpy array can also be used to
        select out members.

        """
        # we can take lists of indices or names; these return a
        # Bundle; repeats already not respected since Bundle functions as a
        # set
        memberlist = self._list()

        if isinstance(index, string_types):
            # a name can be used for indexing
            # a name always returns a Bundle
            out = Bundle([self._list()[i]
                          for i, name in enumerate(self.names)
                          if name == index])

            if not len(out):
                raise KeyError("No name matching string selection")
        elif isinstance(index, slice):
            # we also take slices, obviously
            out = Bundle(*self._list()[index])
            out._cache.update(self._cache)
        else:
            out = super(Bundle, self).__getitem__(index)
        return out

    def __add__(self, other):
        """Addition of collections with collections or treants yields Bundle.

        """
        from .treants import Treant

        if isinstance(other, (Treant, Bundle)):
            return Bundle(self, other)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def __sub__(self, other):
        """Return a Bundle giving the Treants in `a` that are not in `b`.

        Subtracting a Treant from a collection also works.

        """
        from .treants import Treant

        if isinstance(other, Bundle):
            return Bundle(list(set(self) - set(other)))
        elif isinstance(other, Treant):
            return Bundle(list(set(self) - set([other])))
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def __or__(self, other):
        """Return a Bundle giving the union of Bundles `a` and `b`.

        """
        if isinstance(other, Bundle):
            return Bundle(self, other)
        else:
            raise TypeError("Operands must be Bundles.")

    def __and__(self, other):
        """Return a Bundle giving the intersection of Bundles `a` and `b`.

        """
        if isinstance(other, Bundle):
            return Bundle(list(set(self) & set(other)))
        else:
            raise TypeError("Operands must be Bundles.")

    def __xor__(self, other):
        """Return a Bundle giving the symmetric difference of Bundles
        `a` and `b`.

        """
        if isinstance(other, Bundle):
            return Bundle(list(set(self) ^ set(other)))
        else:
            raise TypeError("Operands must be Bundles.")

    def _add(self, *treants):
        """Add any number of members to this collection.

        :Arguments:
            *treants*
                Treants to be added, which may be nested lists/tuples of
                Treants, Bundles, individual Treants or paths to existing
                Treants
        """
        from .treants import Treant

        abspaths = list()
        for treant in treants:
            if treant is None:
                pass
            elif isinstance(treant, (list, tuple, View, Bundle)):
                self._add(*treant)
            elif isinstance(treant, Treant):
                abspaths.append(treant.abspath)
                self._cache[treant.abspath] = treant
            elif isinstance(treant, Tree):
                treantdir = os.path.join(treant.abspath, TREANTDIR_NAME)
                if os.path.exists(treantdir):
                    abspaths.extend(treant.abspath)
                else:
                    raise NotATreantError("Directory '{}' is "
                                          "not a Treant".format(treant))
            elif os.path.exists(treant):
                treantdir = os.path.join(treant, TREANTDIR_NAME)
                if os.path.exists(treantdir):
                    abspaths.append(os.path.abspath(treant))
                else:
                    raise NotATreantError("Directory '{}' is "
                                          "not a Treant".format(treant))
            else:
                raise TypeError("'{}' not a valid input "
                                "for Bundle".format(treant))

        self._add_members(abspaths)

    def _remove(self, *members):
        """Remove any number of members from the collection.

        :Arguments:
            *members*
                instances, indices, names, or absolute paths of the members to
                remove

        """
        from .treants import Treant

        abspaths = self._state
        remove = list()

        for member in members:
            if isinstance(member, int):
                remove.append(abspaths[member])
            elif isinstance(member, Treant):
                remove.append(member.abspath)
            elif isinstance(member, string_types):
                # try abspaths
                abspaths = fnmatch.filter(self.abspaths, member)
                paths = [m.abspath for m in self
                         if m.abspath in abspaths]
                remove.extend(paths)
                # try names
                names = fnmatch.filter(self.names, member)
                paths = [m.abspath for m in self
                         if m.name in names]
                remove.extend(paths)
            else:
                raise TypeError('Only a Treant, index, name, or absolute '
                                'path acceptable')

        self._del_members(remove)

        # remove from cache
        for abspath in remove:
            self._cache.pop(abspath, None)

    @property
    def names(self):
        """Return a list of member names.

        :Returns:
            *names*
                list giving the name of each member, in order

        """
        names = list()
        for member in self._list():
            if member:
                names.append(member.name)
            else:
                names.append(None)

        return names

    @property
    def abspaths(self):
        """Return a list of absolute member directory paths.

        :Returns:
            *abspaths*
                list giving the absolute directory path of each member, in
                order

        """
        return [member.abspath for member in self._list()]

    @property
    def relpaths(self):
        """Return a list of relative member directory paths.

        :Returns:
            *names*
                list giving the relative directory path of each member, in
                order

        """
        return [member.relpath for member in self._list()]

    def _list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the
        collection!

        """
        from .treants import Treant

        abspaths = self._state

        findlist = list()
        memberlist = list()

        for abspath in abspaths:
            if abspath in self._cache:
                memberlist.append(self._cache[abspath])
            elif os.path.exists(os.path.join(abspath, TREANTDIR_NAME)):
                self._cache[abspath] = Treant(abspath)
                memberlist.append(self._cache[abspath])
            else:
                raise NotATreantError("Directory '{}' is "
                                      "not a Treant.".format(abspath))

        return memberlist

    def map(self, function, processes=1, **kwargs):
        """Apply a function to each member, perhaps in parallel.

        A pool of processes is created for *processes* > 1; for example,
        with 40 members and 'processes=4', 4 processes will be created,
        each working on a single member at any given time. When each process
        completes work on a member, it grabs another, until no members remain.

        *kwargs* are passed to the given function when applied to each member

        :Arguments:
            *function*
                function to apply to each member; must take only a single
                treant instance as input, but may take any number of keyword
                arguments

        :Keywords:
            *processes*
                how many processes to use; if 1, applies function to each
                member in member order

        :Returns:
            *results*
                list giving the result of the function for each member,
                in member order; if the function returns ``None`` for each
                member, then only ``None`` is returned instead of a list
            """
        if processes > 1:
            pool = mp.Pool(processes=processes)
            results = dict()
            results = {member.abspath: pool.apply_async(
                    function, args=(member,), kwds=kwargs) for member in self}

            output = {key: results[key].get() for key in results}

            pool.close()
            pool.join()

            # sort by member order
            results = [output[abspath] for abspath in self.abspaths]
        else:
            results = [function(member, **kwargs) for member in self]

        # check if list is all ``None``: if so, we return ``None``
        if all([(i is None) for i in results]):
            results = None

        return results

    def globfilter(self, pattern):
        """Return a Bundle of members that match by name the given globbing
        pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match member names with

        """
        return Bundle([self[name] for name in
                      fnmatch.filter(self.names, pattern)])

    def _add_members(self, abspaths):
        """Add many members at once.

        :Arguments:
            *abspaths*
                list of abspaths

        """
        for abspath in abspaths:
            self._add_member(abspath)

    def _add_member(self, abspath):
        """Add a member to the Bundle.

        :Arguments:
            *abspath*
                absolute path to directory of new member in the filesystem

        """

        if not (abspath in self._state):
            self._state.append(abspath)

    def _del_members(self, abspaths=None, all=False):
        """Remove members from the Bundle.

        :Arguments:
            *abspaths*
                An iterable of abspaths of the members to remove
            *all*
                When True, remove all members [``False``]

        """
        if all:
            self._state = list()
        else:
            for abspath in abspaths:
                try:
                    self._state.remove(abspath)
                    self._cache.pop(abspath, None)
                except ValueError:
                    pass

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, value):
        if isinstance(value, (Tags, list, set)):
            val = list(value)
            self.tags.clear()
            self.tags.add(val)
        else:
            raise TypeError("Can only set with tags, a list, or set")

    @property
    def categories(self):
        return self._categories

    @categories.setter
    def categories(self, value):
        if isinstance(value, (Categories, dict)):
            val = dict(value)
            self.categories.clear()
            self.categories.add(val)
        else:
            raise TypeError("Can only set with categories or dict")

    def get(self, *tags, **categories):
        """Filter to only Treants which match the defined tags and categories.

        If no arguments given, the full Bundle is returned. This method should
        be thought of as a filtering, with more values specified giving only
        those Treants that match.

        Parameters
        ----------
        *tags
            Tags to match.
        **categories
            Category key, value pairs to match.

        Returns
        -------
        Bundle
            All matched Treants.

        Examples
        --------
        Doing a `get` with::

        >>> b.get('this')  # doctest: +SKIP

        is equivalent to::

        >>> b.tags.filter('this')  # doctest: +SKIP

        Finally, doing::

        >>> b.get('this', length=5)  # doctest: +SKIP

        is equivalent to::

        >>> b_n = b.tags.filter('this')  # doctest: +SKIP
        >>> b_n.categories.groupby('length')[5.0]  # doctest: +SKIP

        """
        if not (tags or categories):
            # if nothing given, return an empty version of this bundle
            return self

        output = self  # initially return self
        for c, v in categories.items():
            try:
                output &= output.categories.groupby(c)[v]
            except KeyError:
                # the key wasn't present, so we return an empty Bundle
                output = Bundle()
                break  # don't need to loop any further

        # if any tags are given that don't match anything, this will yield an
        # empty Bundle
        for t in tags:
            output &= output.tags.filter(t)

        return output


class _Loc(object):
    """Path accessor for collections."""

    def __init__(self, collection):
        self._collection = collection

    def __getitem__(self, path):
        """Get Tree/Leaf at `path` relative to each Tree in collection.

        """
        return View([t[path] for t in self._collection
                     if isinstance(t, Tree)])


class _TreeLoc(_Loc):
    """Tree accessor for collections."""

    def __getitem__(self, path):
        """Get Tree at `path` relative to each Tree in collection.

        """
        return View([t.treeloc[path] for t in self._collection
                     if isinstance(t, Tree)])


class _LeafLoc(_Loc):
    """Leaf accessor for collections."""

    def __getitem__(self, path):
        """Get Leaf at `path` relative to each Tree in collection.

        """
        return View([t.leafloc[path] for t in self._collection
                     if isinstance(t, Tree)])
