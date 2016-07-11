"""
The Bundle object is the primary manipulator for Treants in aggregate. They are
returned as queries to Groups and other Bundles. They offer convenience methods
for dealing with many Treants at once. Views give the same kind of aggregation
conveniences for Trees and Leaves.

"""

from __future__ import absolute_import

import os
import functools
from collections import namedtuple, defaultdict

import multiprocessing as mp
import glob
import fnmatch

from six import string_types
from six.moves import zip

from . import filesystem
from . import _AGGLIMBS, _AGGTREELIMBS
from .trees import Tree, Leaf
from .manipulators import discover


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

    def glob(self, pattern):
        """Return a View of all child Leaves and Trees of members matching
        given globbing pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match files and directories with

        """
        return View([member.glob(pattern) for member in self
                     if isinstance(member, Tree)], limbs=self.limbs)

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
    def limbs(self):
        """A set giving the names of this collection's attached limbs.

        """
        return self._classagglimbs | self._agglimbs


class View(CollectionMixin):
    """An ordered set of Trees and Leaves.

    Parameters
    ----------
    vegs : Tree, Leaf, or list
        Trees and/or Leaves to be added, which may be nested lists of Trees
        and Leaves. Trees and Leaves can be given as either objects or
        paths.
    limbs : list or set
        Names of limbs to immediately attach.

    """
    _classagglimbs = set()
    _agglimbs = set()

    def __init__(self, *vegs, **kwargs):
        self._state = list()
        self.add(*vegs)

        # attach any limbs given
        for agglimb in kwargs.pop('limbs', []):
            try:
                self.attach(agglimb)
            except KeyError:
                pass

    def __repr__(self):
        return "<View({})>".format(self._list())

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        A single integer will yield a single member. Lists of either will yield
        a View with members in order by the given items. Giving a basename will
        always yield a View, since more than one Tree or Leaf may have the same
        basename and so they are not guaranteed to be unique.

        A boolean index by way of a list or numpy array can also be used to
        select out members.

        """
        # we can take lists of indices, names, or uuids; these return a
        # View; repeats already not respected since View functions as a
        # set
        if ((isinstance(index, list) or hasattr(index, 'dtype')) and
                all([isinstance(item, bool) for item in index])):
            # boolean indexing
            memberlist = self._list()
            out = View([memberlist[i] for i, val in enumerate(index) if val],
                       limbs=self.limbs)
        elif isinstance(index, list):
            memberlist = self._list()
            out = View([memberlist[item] for item in index],
                       limbs=self.limbs)
        elif isinstance(index, int):
            # an index gets the member at that position
            out = self._list()[index]
        elif isinstance(index, string_types):
            # a name can be used for indexing
            # always returns a View
            out = View([self._list()[i] for i, name
                        in enumerate(self.names) if name == index],
                       limbs=self.limbs)

        elif isinstance(index, slice):
            # we also take slices, obviously
            out = View(*self._list()[index], limbs=self.limbs)
        else:
            raise IndexError("Cannot index View with given values")

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
        if isinstance(other, (Tree, Leaf, View, list)):
            limbs = self.limbs | other.limbs
            return View(self, other, limbs=limbs)
        else:
            raise TypeError("Right operand must be a Tree, Leaf, or View.")

    def __sub__(self, other):
        """Return a View giving the Treants in `self` that are not in `other`.

        Subtracting a Tree/Leaf from a View also works.

        """
        if isinstance(other, View):
            limbs = self.limbs | other.limbs
            return View(list(set(self) - set(other)), limbs=limbs)
        elif isinstance(other, (Tree, Leaf)):
            limbs = self.limbs | other.limbs
            return View(list(set(self) - set([other])), limbs=limbs)
        else:
            raise TypeError("Right operand must be a Tree, Leaf, or View.")

    def __or__(self, other):
        """Return a View giving the union of Views `self` and `other`.

        """
        if isinstance(other, View):
            limbs = self.limbs | other.limbs
            return View(self, other, limbs=limbs)
        else:
            raise TypeError("Operands must be Views.")

    def __and__(self, other):
        """Return a View giving the intersection of Views `self` and `other`.

        """
        if isinstance(other, View):
            limbs = self.limbs | other.limbs
            return View(list(set(self) & set(other)), limbs=limbs)
        else:
            raise TypeError("Operands must be Views.")

    def __xor__(self, other):
        """Return a View giving the symmetric difference of Views `self` and
        `other`.

        """
        if isinstance(other, View):
            limbs = self.limbs | other.limbs
            return View(list(set(self) ^ set(other)), limbs=limbs)
        else:
            raise TypeError("Operands must be Views.")

    @classmethod
    def _attach_aggtreelimb_class(cls, limb):
        """Attach an aggtreelimb to the class.

        """
        # property definition
        def getter(self):
            if not hasattr(self, "_"+limb._name):
                setattr(self, "_"+limb._name, limb(self))
            return getattr(self, "_"+limb._name)

        # set the property
        setattr(cls, limb._name,
                property(getter, None, None, limb.__doc__))

        if limb._name in _AGGTREELIMBS:
            cls._classagglimbs.add(limb._name)

    def _attach_aggtreelimb(self, limb):
        """Attach an aggtreelimb.

        """
        try:
            setattr(self, limb._name, limb(self))
        except AttributeError:
            pass

        if limb._name in _AGGTREELIMBS:
            self._agglimbs.add(limb._name)

    def attach(self, *aggtreelimbname):
        """Attach aggtreelimbs by name to this View. Attaches corresponding limb
        to any member Trees.

        """
        for ln in aggtreelimbname:
            # try and get the aggtreelimb class specified
            try:
                aggtreelimb = _AGGTREELIMBS[ln]
            except KeyError:
                raise KeyError("No such aggtreelimb '{}'".format(ln))

            # attach agglimb; if it's already there, that's okay
            try:
                self._attach_aggtreelimb(aggtreelimb)
            except AttributeError:
                pass

            # attach limb to each member
            for member in self._list():
                if isinstance(member, Tree):
                    member.attach(ln)

    def add(self, *vegs):
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
                self.add(*veg)
            elif isinstance(veg, (View, Bundle)):
                self.add(*list(veg))
            elif isinstance(veg, Treant):
                outconts.append(veg.tree)
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
        return View([member for member in self if isinstance(member, Tree)],
                    limbs=self.limbs)

    @property
    def memberleaves(self):
        """A View giving only members that are Leaves (or subclasses).

        """
        return View([member for member in self if isinstance(member, Leaf)],
                    limbs=self.limbs)

    @property
    def children(self):
        """A View of all children within the member Trees.

        """
        return View([member.children for member in self.membertrees],
                    limbs=self.limbs)

    @property
    def trees(self):
        """A View of directories within the member Trees.

        Hidden directories are not included.

        """
        return View([member.trees for member in self.membertrees],
                    limbs=self.limbs)

    @property
    def leaves(self):
        """A View of the files within the member Trees.

        Hidden files are not included.

        """
        return View([member.leaves for member in self.membertrees],
                    limbs=self.limbs)

    @property
    def hidden(self):
        """A View of the hidden files and directories within the member Trees.

        """
        return View([member.hidden for member in self.membertrees],
                    limbs=self.limbs)

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
    def bundle(self):
        """Obtain a Bundle of all existing Treants among the Trees and Leaves
        in this View.

        """
        return Bundle(self, limbs=self.limbs)

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
        """Return a View of members that match by name the given globbing
        pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match member names with

        """
        return View([self[name] for name in
                     fnmatch.filter(self.names, pattern)], limbs=self.limbs)

    def make(self):
        """Make the Trees and Leaves in this View if they don't already exist.

        Returns
        -------
        View
            This View.

        """
        for member in self:
            self.make()

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
    limbs : list or set
        Names of limbs to immediately attach.

    """
    _memberpaths = ['abspath']
    _fields = ['uuid', 'treanttype']
    _fields.extend(_memberpaths)

    _classagglimbs = set()
    _agglimbs = set()

    def __init__(self, *treants, **kwargs):
        self._cache = dict()
        self._state = list()
        self._searchtime = 10

        self.add(*treants)

        # attach any limbs given
        for agglimb in kwargs.pop('limbs', []):
            if agglimb not in self.limbs:
                try:
                    self.attach(agglimb)
                except KeyError:
                    pass

    def __repr__(self):
        return "<Bundle({})>".format(self._list())

    def __str__(self):
        out = "<- Bundle ->\n"

        for member in self._list():
            out += "  {}\n".format(member)

        out += "<- ---- ->"
        return out

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        A single integer or uuid will yield a single Treant. Lists of
        either will yield a Bundle with members in order by the given items.
        Giving a name will always yield a Bundle, since names are not
        guaranteed to be unique.

        A boolean index by way of a list or numpy array can also be used to
        select out members.

        """
        # we can take lists of indices, names, or uuids; these return a
        # Bundle; repeats already not respected since Bundle functions as a
        # set
        if ((isinstance(index, list) or hasattr(index, 'dtype')) and
                all([isinstance(item, bool) for item in index])):
            # boolean indexing
            memberlist = self._list()
            out = Bundle([memberlist[i] for i, val in enumerate(index) if val],
                         limbs=self.limbs)
        elif isinstance(index, list):
            memberlist = self._list()
            out = Bundle([memberlist[item] for item in index],
                         limbs=self.limbs)
        elif isinstance(index, int):
            # an index gets the member at that position
            out = self._list()[index]
        elif isinstance(index, string_types):
            # a name or uuid can be used for indexing
            # a name always returns a Bundle
            out = Bundle([self.filepaths[i] for i, name
                          in enumerate(self.names) if name == index],
                         limbs=self.limbs)

            # if no names match, we try uuids
            if not len(out):
                out = [member for member in self if member.uuid == index]
                if not len(out):
                    raise KeyError("No name or uuid matching string selection")
                else:
                    # we want to return a Treant, not a list for uuid matches
                    out = out[0]
        elif isinstance(index, slice):
            # we also take slices, obviously
            out = Bundle(*self.filepaths[index], limbs=self.limbs)
            out._cache.update(self._cache)
        else:
            raise IndexError("Cannot index Bundle with given values")

        return out

    def __add__(self, other):
        """Addition of collections with collections or treants yields Bundle.

        """
        from .treants import Treant

        if isinstance(other, (Treant, Bundle, list)):
            limbs = self.limbs | other.limbs
            return Bundle(self, other, limbs=limbs)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def __sub__(self, other):
        """Return a Bundle giving the Treants in `a` that are not in `b`.

        Subtracting a Treant from a collection also works.

        """
        from .treants import Treant

        if isinstance(other, Bundle):
            limbs = self.limbs | other.limbs
            return Bundle(list(set(self) - set(other)), limbs=limbs)
        elif isinstance(other, Treant):
            limbs = self.limbs | other.limbs
            return Bundle(list(set(self) - set([other])), limbs=limbs)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def __or__(self, other):
        """Return a Bundle giving the union of Bundles `a` and `b`.

        """
        if isinstance(other, Bundle):
            limbs = self.limbs | other.limbs
            return Bundle(self, other, limbs=limbs)
        else:
            raise TypeError("Operands must be Bundles.")

    def __and__(self, other):
        """Return a Bundle giving the intersection of Bundles `a` and `b`.

        """
        if isinstance(other, Bundle):
            limbs = self.limbs | other.limbs
            return Bundle(list(set(self) & set(other)), limbs=limbs)
        else:
            raise TypeError("Operands must be Bundles.")

    def __xor__(self, other):
        """Return a Bundle giving the symmetric difference of Bundles
        `a` and `b`.

        """
        if isinstance(other, Bundle):
            limbs = self.limbs | other.limbs
            return Bundle(list(set(self) ^ set(other)), limbs=limbs)
        else:
            raise TypeError("Operands must be Bundles.")

    @classmethod
    def _attach_agglimb_class(cls, limb):
        """Attach a agglimb to the class.

        """
        # property definition
        def getter(self):
            if not hasattr(self, "_"+limb._name):
                setattr(self, "_"+limb._name, limb(self))
            return getattr(self, "_"+limb._name)

        try:
            setter = limb._setter
        except AttributeError:
            setter = None

        # set the property
        setattr(cls, limb._name,
                property(getter, setter, None, limb.__doc__))

        if limb._name in _AGGTREELIMBS or limb._name in _AGGLIMBS:
            cls._classagglimbs.add(limb._name)

    def _attach_agglimb(self, limb):
        """Attach an agglimb.

        """
        try:
            setattr(self, limb._name, limb(self))
        except AttributeError:
            pass

        if limb._name in _AGGTREELIMBS or limb._name in _AGGLIMBS:
            self._agglimbs.add(limb._name)

    def attach(self, *agglimbname):
        """Attach agglimbs by name to this collection. Attaches corresponding limb
        to member Treants.

        """
        for ln in agglimbname:
            # try and get the aggtreelimb class specified
            try:
                agglimb = _AGGTREELIMBS[ln]
            except KeyError:
                # if not an aggtreelimb, perhaps its an agglimb?
                try:
                    agglimb = _AGGLIMBS[ln]
                except KeyError:
                    raise KeyError("No such agglimb '{}'".format(ln))

            # attach agglimb; if it's already there, that's okay
            try:
                self._attach_agglimb(agglimb)
            except AttributeError:
                pass

            # attach limb to each member
            for member in self._list():
                member.attach(ln)

    def add(self, *treants):
        """Add any number of members to this collection.

        :Arguments:
            *treants*
                treants to be added, which may be nested lists of treants;
                treants can be given as either objects or paths to directories
                that contain treant statefiles; glob patterns are also allowed,
                and all found treants will be added to the collection
        """
        from .treants import Treant

        outconts = list()
        for treant in treants:
            if treant is None:
                pass
            elif isinstance(treant, (list, tuple, View)):
                self.add(*treant)
            elif isinstance(treant, Bundle):
                self.add(*treant.filepaths)
                self._cache.update(treant._cache)
            elif isinstance(treant, Treant):
                outconts.append(treant)
                self._cache[treant.uuid] = treant
            elif isinstance(treant, (Leaf, Tree)):
                tre = filesystem.path2treant(treant.abspath)
                outconts.extend(tre)
            elif os.path.exists(treant):
                tre = filesystem.path2treant(treant)
                outconts.extend(tre)
            elif isinstance(treant, string_types):
                tre = filesystem.path2treant(*glob.glob(treant))
                outconts.extend(tre)
            else:
                raise TypeError("'{}' not a valid input "
                                "for Bundle".format(treant))

        attrs = []
        for attr in ('uuid', 'treanttype', 'abspath'):
            attrs.append([getattr(treant, attr) for treant in outconts])

        self._add_members(*attrs)

    def remove(self, *members):
        """Remove any number of members from the collection.

        :Arguments:
            *members*
                instances or indices of the members to remove

        """
        from .treants import Treant

        uuids = self._get_members_uuid()
        remove = list()
        for member in members:
            if isinstance(member, int):
                remove.append(uuids[member])
            elif isinstance(member, Treant):
                remove.append(member.uuid)
            elif isinstance(member, string_types):
                names = fnmatch.filter(self.names, member)
                uuids = [member.uuid for member in self
                         if member.name in names]
                remove.extend(uuids)

            else:
                raise TypeError('Only an integer or treant acceptable')

        self._del_members(remove)

        # remove from cache
        for uuid in remove:
            self._cache.pop(uuid, None)

    def clear(self):
        """Remove all members.

        """
        self._del_members(all=True)

    @property
    def treanttypes(self):
        """Return a list of member treanttypes.

        """
        return self._get_members_treanttype()

    @property
    def names(self):
        """Return a list of member names.

        Members that can't be found will have name ``None``.

        :Returns:
            *names*
                list giving the name of each member, in order;
                members that are missing will have name ``None``

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

        Members that can't be found will have path ``None``.

        :Returns:
            *names*
                list giving the absolute directory path of each member, in
                order; members that are missing will have path ``None``

        """
        return [member.abspath if member else None for member in self._list()]

    @property
    def relpaths(self):
        """Return a list of relative member directory paths.

        Members that can't be found will have path ``None``.

        :Returns:
            *names*
                list giving the relative directory path of each member, in
                order; members that are missing will have path ``None``

        """
        return [member.relpath if member else None for member in self._list()]

    @property
    def filepaths(self):
        """Return a list of member filepaths.

        Members that can't be found will have filepath ``None``.

        :Returns:
            *names*
                list giving the filepath of each member, in order;
                members that are missing will have filepath ``None``

        """
        filepaths = list()
        for member in self._list():
            if member:
                filepaths.append(member.filepath)
            else:
                filepaths.append(None)

        return filepaths

    @property
    def uuids(self):
        """Return a list of member uuids.

        :Returns:
            *uuids*
                list giving the uuid of each member, in order

        """
        return self._get_members_uuid()

    def _check(self):
        """Check that all member paths resolve properly.

        """
        members = self._get_members()

        paths = {path: members[path] for path in self._memberpaths}

        foxhound = filesystem.Foxhound(self, members['uuid'], paths,
                                       timeout=self.searchtime)
        found = foxhound.fetch(as_treants=False)

        if None not in found.values():
            treanttypes = [os.path.basename(path).split(os.extsep)[0]
                           for path in found.values()]

            self._add_members(found.keys(), treanttypes, found.values())
            return True
        else:
            return False

    def _list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the
        collection!

        Missing members will be present in the list as ``None``. This method is
        not intended for user-level use.

        """
        members = self._get_members()
        uuids = members['uuid']

        findlist = list()
        memberlist = list()

        for uuid in uuids:
            if uuid in self._cache and self._cache[uuid]:
                memberlist.append(self._cache[uuid])
            else:
                memberlist.append(None)
                findlist.append(uuid)

        # track down our non-cached treants
        if findlist:
            paths = {path: members[path]
                     for path in self._memberpaths}
            foxhound = filesystem.Foxhound(self, findlist, paths,
                                           timeout=self.searchtime)
            foundconts = foxhound.fetch(as_treants=True)

            # add to cache, and ensure we get updated paths with a re-add in
            # case of an IOError, skip (probably due to permissions, but will
            # need something more robust later
            self._cache.update(foundconts)
            try:
                self.add(*foundconts.values())
            except OSError:
                pass

            # insert found treants into output list
            for uuid in findlist:
                result = foundconts[uuid]
                if not result:
                    ind = list(members['uuid']).index(uuid)
                    raise IOError("Could not find member {} (uuid: {});"
                                  " re-add or remove it.".format(ind, uuid))

                memberlist[list(uuids).index(uuid)] = result

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
            results = {member.uuid: pool.apply_async(
                    function, args=(member,), kwds=kwargs) for member in self}

            output = {key: results[key].get() for key in results}

            pool.close()
            pool.join()

            # sort by member order
            results = [output[uuid] for uuid in self.uuids]
        else:
            results = [function(member, **kwargs) for member in self]

        # check if list is all ``None``: if so, we return ``None``
        if all([(i is None) for i in results]):
            results = None

        return results

    @property
    def searchtime(self):
        """Max time to spend searching for missing members, in seconds.

        Setting a larger value allows more time for the collection to look for
        members elsewhere in the filesystem.

        If `None`, there will be no time limit. Use with care.

        """
        return self._searchtime

    @searchtime.setter
    def searchtime(self, value):
        if isinstance(value, (float, int)) or value is None:
            self._searchtime = value
        else:
            raise TypeError("Must give a number or `None` for searchtime")

    def flatten(self, exclude=None):
        """Return a flattened version of this Bundle.

        The resulting Bundle will have all members of any member Groups,
        without the Groups.

        Parameters
        ----------
        exclude : list
            uuids of Groups to leave out of flattening; these will not in the
            resulting Bundle.

        Returns
        -------
        flattened : Bundle
            the flattened Bundle with no Groups

        """
        if not exclude:
            exclude = list()

        guuids = list(exclude)
        memberlist = self._list()
        flattened = Bundle(limbs=self.limbs)

        for member in memberlist:
            if hasattr(member, 'members') and member.uuid not in exclude:
                guuids.append(member.uuid)
                flattened += member.members.flatten(guuids)
            elif not hasattr(member, 'members'):
                flattened.add(member)

        return flattened

    @property
    def view(self):
        """Obtain a View giving the Tree for each Treant in this Bundle.

        """
        return View([member.tree for member in self], limbs=self.limbs)

    def globfilter(self, pattern):
        """Return a Bundle of members that match by name the given globbing
        pattern.

        Parameters
        ----------
        pattern : string
            globbing pattern to match member names with

        """
        return Bundle([self[name] for name in
                      fnmatch.filter(self.names, pattern)], limbs=self.limbs)

    def _add_members(self, uuids, treanttypes, abspaths):
        """Add many members at once.

        Given lists must be in the same order with respect to the members they
        describe.

        :Arguments:
            *uuids*
                list of uuids
            *treanttypes*
                list of treanttypes
            *abspaths*
                list of abspaths

        """
        for uuid, treanttype, abspath in zip(uuids, treanttypes, abspaths):
            self._add_member(uuid, treanttype, abspath)

    def _add_member(self, uuid, treanttype, abspath):
        """Add a member to the Bundle.

        If the member is already present, its location will be updated with
        the given location.

        :Arguments:
            *uuid*
                the uuid of the new member
            *treanttype*
                the treant type of the new member
            *abspath*
                absolute path to directory of new member in the filesystem

        """
        member_rec = {'uuid': uuid,
                      'treanttype': treanttype,
                      'abspath': os.path.abspath(abspath)}

        # check if uuid already present
        uuids = [member['uuid'] for member in self._state]

        if uuid in uuids:
            self._state[uuids.index(uuid)] = member_rec
        else:
            self._state.append(member_rec)

    def _del_members(self, uuids=None, all=False):
        """Remove members from the Bundle.

        :Arguments:
            *uuids*
                An iterable of uuids of the members to remove
            *all*
                When True, remove all members [``False``]

        """
        if all:
            self._state = list()
        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uuid) for uuid in uuids])

            # get matching rows
            memberlist = list()
            for i, member in enumerate(self._state):
                for uuid in uuids:
                    if (member['uuid'] == uuid):
                        memberlist.append(i)

            memberlist.sort()
            j = 0
            # delete matching entries; have to use j to shift the register as
            # we remove entries
            for i in memberlist:
                self._state.pop(i - j)
                j = j + 1

    def _get_member(self, uuid):
        """Get all stored information on the specified member.

        Returns a dictionary whose keys are column names and values the
        corresponding values for the member.

        :Arguments:
            *uuid*
                uuid of the member to retrieve information for

        :Returns:
            *memberinfo*
                a dictionary containing all information stored for the
                specified member
        """
        memberinfo = None
        for member in self._state:
            if member['uuid'] == uuid:
                memberinfo = member

        return memberinfo

    def _get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                dict giving full member data, with fields as keys and in member
                order
        """
        out = defaultdict(list)

        for member in self._state:
            for key in self._fields:
                out[key].append(member[key])

        return out

    def _get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving treanttype of each member, in order
        """
        return [member['uuid'] for member in self._state]

    def _get_members_names(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving treanttype of each member, in order
        """
        return [os.path.basename(member['abspath']) for member in self._state]

    def _get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        return [member['treanttype'] for member in self._state]


class _Loc(object):
    """Subtree accessor for collections."""

    def __init__(self, collection):
        self._collection = collection

    def __getitem__(self, path):
        """Get Tree/Leaf at `path` relative to each Tree in collection.

        """
        return View([t[path] for t in self._collection
                     if isinstance(t, Tree)], limbs=self._collection.limbs)
