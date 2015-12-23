"""
The Bundle object is the primary manipulator for Treants in aggregate.
They are returned as queries to Groups, Coordinators, and other Bundles. They
offer convenience methods for dealing with many Treants at once.

"""
import os

import multiprocessing as mp
import glob
import fnmatch

from datreant import backends
from datreant import filesystem
import datreant.treants


class CollectionBase(object):
    """Common interface elements for ordered sets of Treants.

    :class:`datreant.limbs.Members` and :class:`Bundle` both use this
    interface.

    """
    def __len__(self):
        return len(self._list())

    @classmethod
    def _attach_limb(cls, limb):
        """Attach a limb to the class, or to self if an instance.

        """
        # property definition
        def getter(self):
            if not hasattr(self, "_"+limb._name):
                setattr(self, "_"+limb._name, limb(self))
            return getattr(self, "_"+limb._name)

        # set the property
        setattr(cls, limb._name,
                property(getter, None, None, limb.__doc__))

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        """
        if isinstance(index, int):
            out = self._list()[index]
        else:
            out = Bundle(*self._list()[index])

        return out

    def __add__(a, b):
        """Addition of collections with collections or treants yields Bundle.

        """
        if (isinstance(a, (datreant.treants.Treant, CollectionBase)) and
                isinstance(b, (datreant.treants.Treant, CollectionBase))):
            return Bundle(a, b)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def add(self, *treants):
        """Add any number of members to this collection.

        :Arguments:
            *treants*
                treants to be added, which may be nested lists of treants;
                treants can be given as either objects or paths to directories
                that contain treant statefiles; glob patterns are also allowed,
                and all found treants will be added to the collection
        """
        from datreant.limbs import Members
        from datreant.treants import Treant

        outconts = list()
        for treant in treants:
            if treant is None:
                pass
            elif isinstance(treant,
                            (list, tuple, CollectionBase)):
                self.add(*treant)
            elif isinstance(treant, Treant):
                outconts.append(treant)
            elif os.path.exists(treant):
                tre = filesystem.path2treant(treant)
                for t in tre:
                    outconts.append(t)
            else:
                tre = filesystem.path2treant(*glob.glob(treant))
                for t in tre:
                    outconts.append(t)

        for treant in outconts:
            self._backend.add_member(treant.uuid,
                                     treant.treanttype,
                                     treant.basedir)

    def remove(self, *members, **kwargs):
        """Remove any number of members from the Group.

        :Arguments:
            *members*
                instances or indices of the members to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        from .treants import Treant

        uuids = self._backend.get_members_uuid()
        if kwargs.pop('all', False):
            remove = uuids
        else:
            remove = list()
            for member in members:
                if isinstance(member, int):
                    remove.append(uuids[member])
                elif isinstance(member, Treant):
                    remove.append(member.uuid)
                elif isinstance(member, basestring):
                    names = fnmatch.filter(self.names, member)
                    uuids = [member.uuid for member in self
                             if (member.name in names)]
                    remove.extend(uuids)

                else:
                    raise TypeError('Only an integer or treant acceptable')

        self._backend.del_member(*remove)

        # remove from cache
        for uuid in remove:
            self._cache.pop(uuid, None)

    @property
    def treanttypes(self):
        """Return a list of member treanttypes.

        """
        return self._backend.get_members_treanttype()

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
    def basedirs(self):
        """Return a list of member basedirs.

        Members that can't be found will have basedir ``None``.

        :Returns:
            *names*
                list giving the basedir of each member, in order;
                members that are missing will have basedir ``None``

        """
        basedirs = list()
        for member in self._list():
            if member:
                basedirs.append(member.basedir)
            else:
                basedirs.append(None)

        return basedirs

    @property
    def filepath(self):
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
        return self._backend.get_members_uuid()

    def _list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        Missing members will be present in the list as ``None``. This method is
        not intended for user-level use.

        """
        members = self._backend.get_members()
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
        paths = {path: members[path]
                 for path in self._backend.memberpaths}
        foxhound = filesystem.Foxhound(self, findlist, paths)
        foundconts = foxhound.fetch(as_treants=True)

        # add to cache, and ensure we get updated paths with a re-add
        # in case of an IOError, skip (probably due to permissions, but will
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
                raise IOError("Could not find member" +
                              " {} (uuid: {});".format(ind, uuid) +
                              " re-add or remove it.")

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
            for member in self:
                results[member.uuid] = pool.apply_async(
                        function, args=(member,), kwds=kwargs).get()
            pool.close()
            pool.join()

            # sort by member order
            results = [results[uuid] for uuid in self.uuids]
        else:
            results = [function(member, **kwargs) for member in self]

        # check if list is all ``None``: if so, we return ``None``
        if all([(i is None) for i in results]):
            results = None

        return results


class _BundleBackend():
    """Backend class for Bundle.

    Has same interface as Group-specific components of
    :class:`backends.GroupFile`. Behaves practically like an in-memory
    version of a state-file, but with only the components needed for the
    Bundle.

    """
    memberpaths = ['abspath']
    fields = ['uuid', 'treanttype', 'abspath']

    def __init__(self):
        self.record = list()

    def add_member(self, uuid, treanttype, basedir):
        """Add a member to the Bundle.

        If the member is already present, its location will be updated with
        the given location.

        :Arguments:
            *uuid*
                the uuid of the new member
            *treanttype*
                the treant type of the new member
            *basedir*
                basedir of the new member in the filesystem

        """
        # check if uuid already present
        uuids = [member[0] for member in self.record]

        if uuid not in uuids:
            self.record.append([uuid,
                                treanttype,
                                os.path.abspath(basedir)])

    def del_member(self, *uuid, **kwargs):
        """Remove a member from the Group.

        :Arguments:
            *uuid*
                the uuid(s) of the member(s) to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        purge = kwargs.pop('all', False)

        if purge:
            self.record = list()
        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uid) for uid in uuid])

            # get matching rows
            # TODO: possibly faster to use table.where
            memberlist = list()
            for i, member in enumerate(self.record):
                for uuid in uuids:
                    if (member[0] == uuid):
                        memberlist.append(i)

            memberlist.sort()
            j = 0
            # delete matching entries; have to use j to shift the register as
            # we remove entries
            for i in memberlist:
                self.record.pop(i - j)
                j = j + 1

    def get_member(self, uuid):
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
        for member in self.record:
            if member[0] == uuid:
                memberinfo = member

        if memberinfo:
            memberinfo = {x: y for x, y in zip(self.fields, memberinfo)}

        return memberinfo

    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                dict giving full member data, with fields as keys and in member
                order
        """
        out = {key: [] for key in self.fields}

        for member in self.record:
            for i, key in enumerate(self.fields):
                out[key].append(member[i])

        return out

    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving treanttype of each member, in order
        """
        return [member[0] for member in self.record]

    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                list giving treanttype of each member, in order
        """
        return [member[1] for member in self.record]

    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                list containing all paths to member basedirs, in member order
        """
        return [member[2:] for member in self.record]


class Bundle(CollectionBase):
    """Non-persistent collection of treants.

    A Bundle is basically an indexable set. It is often used to return the
    results of a query on a Coordinator or a Group, but can be used on its
    own as well.

    """

    def __init__(self, *treants, **kwargs):
        """Generate a Bundle from any number of Treants.

        :Arguments:
            *treants*
                treants to be added, which may be nested lists of treants;
                treants can be given as either objects or paths to directories
                that contain treant statefiles; glob patterns are also allowed,
                and all found treants will be added to the collection
        """
        self._backend = _BundleBackend()
        self._cache = dict()

        self.add(*treants)

    def __repr__(self):
        return "<Bundle({})>".format(self._list())
