"""
The Bundle object is the primary manipulator for Treants in aggregate.
They are returned as queries to Groups, Coordinators, and other Bundles. They
offer convenience methods for dealing with many Treants at once.

"""
import os

import numpy as np
import multiprocessing as mp

from datreant import persistence
from datreant import filesystem


class _CollectionBase(object):
    """Common interface elements for ordered sets of Treants.

    :class:`aggregators.Members` and :class:`Bundle` both use this interface.

    """
    def __len__(self):
        return len(self._list())

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.

        """
        if isinstance(index, int):
            out = self._list()[index]
        else:
            out = Bundle(*self._list()[index])

        return out

    def add(self, *treants):
        """Add any number of members to this collection.

        :Arguments:
            *treants*
                Treants and/or Groups to be added; may be a list of Treants
                and/or Groups; Treants or Groups can be given as either objects
                or paths to directories that contain object statefiles
        """
        from datreant.aggregators import Members
        from datreant.treants import Treant

        outconts = list()
        for treant in treants:
            if treant is None:
                pass
            elif isinstance(treant,
                            (list, tuple, Bundle, Members)):
                self.add(*treant)
            elif isinstance(treant, Treant):
                outconts.append(treant)
            elif os.path.exists(treant):
                cont = filesystem.path2treant(treant)
                for c in cont:
                    outconts.append(c)

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
        return self._backend.get_members_treanttype().tolist()

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
    def uuids(self):
        """Return a list of member uuids.

        :Returns:
            *uuids*
                list giving the uuid of each member, in order

        """
        return self._backend.get_members_uuid().tolist()

    def _list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        Missing members will be present in the list as ``None``. This method is
        not intended for user-level use.

        """
        members = self._backend.get_members()
        uuids = members['uuid'].flatten().tolist()

        findlist = list()
        memberlist = list()

        for uuid in uuids:
            if uuid in self._cache and self._cache[uuid]:
                memberlist.append(self._cache[uuid])
            else:
                memberlist.append(None)
                findlist.append(uuid)

        # track down our non-cached treants
        paths = {path: members[path].flatten().tolist()
                 for path in self._backend.memberpaths}
        foxhound = filesystem.Foxhound(self, findlist, paths)
        foundconts = foxhound.fetch(as_treants=True)

        # add to cache, and ensure we get updated paths with a re-add
        # in case of an IOError, skip (probably due to permissions, but will
        # need something more robust later
        self._cache.update(foundconts)
        try:
            self.add(*foundconts.values())
        except IOError:
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

    @property
    def data(self):
        """Access the data of each member, collectively.

        """
        from .aggregators import MemberData
        if not self._data:
            self._data = MemberData(self)
        return self._data

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
    :class:`persistence.GroupFile`. Behaves practically like an in-memory
    version of a state-file, but with only the components needed for the
    Bundle.

    """
    memberpaths = ['abspath']

    def __init__(self):
        # our table will be a structured array matching the schema of the
        # GroupFile _Members Table
        self.table = np.array(
                [],
                dtype={'names': ['uuid', 'treanttype', 'abspath'],
                       'formats': ['a{}'.format(persistence.uuidlength),
                                   'a{}'.format(persistence.namelength),
                                   'a{}'.format(persistence.pathlength)]
                       }).reshape(1, -1)

    def _member2record(self, uuid, treanttype, basedir):
        """Return a record array from a member's information.

        This method defines the scheme for the Bundle's record array.

        """
        return np.array(
                (uuid, treanttype, os.path.abspath(basedir)),
                dtype={'names': ['uuid', 'treanttype', 'abspath'],
                       'formats': ['a{}'.format(persistence.uuidlength),
                                   'a{}'.format(persistence.namelength),
                                   'a{}'.format(persistence.pathlength)]
                       }).reshape(1, -1)

    def add_member(self, uuid, treanttype, basedir):
        """Add a member to the Group.

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
        if self.table.shape == (1, 0):
            self.table = self._member2record(uuid, treanttype, basedir)
        else:
            # check if uuid already present
            index = np.where(self.table['uuid'] == uuid)[0]
            if index.size > 0:
                # if present, update location
                self.table[index[0]]['abspath'] = os.path.abspath(basedir)
            else:
                newmem = self._member2record(uuid, treanttype, basedir)
                self.table = np.vstack((self.table, newmem))

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
            self.__init__()
        else:
            # remove redundant uuids from given list if present
            uuids = set([str(uid) for uid in uuid])

            # remove matching elements
            matches = list()
            for uuid in uuids:
                index = np.where(self.table['uuid'] == uuid)[0]
                if index:
                    matches.append(index)

            self.table = np.delete(self.table, matches)

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
        memberinfo = self.table[self.table[uuid] == uuid]

        if memberinfo:
            memberinfo = {x: memberinfo[x] for x in memberinfo.dtype.names}
        else:
            memberinfo = None

        return memberinfo

    def get_members(self):
        """Get full member table.

        Sometimes it is useful to read the whole member table in one go instead
        of doing multiple reads.

        :Returns:
            *memberdata*
                structured array giving full member data, with
                each row corresponding to a member
        """
        return self.table

    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                array giving treanttype of each member, in order
        """
        return self.table['uuid'].flatten()

    def get_members_treanttype(self):
        """List treanttype for each member.

        :Returns:
            *treanttypes*
                array giving treanttype of each member, in order
        """
        return self.table['treanttype'].flatten()

    def get_members_basedir(self):
        """List basedir for each member.

        :Returns:
            *basedirs*
                structured array containing all paths to member basedirs
        """
        return self.table['abspath'].flatten()


class Bundle(_CollectionBase):
    """Non-persistent Treant for Treants and Groups.

    A Bundle is basically an indexable set. It is often used to return the
    results of a query on a Coordinator or a Group, but can be used on its
    own as well.

    """

    def __init__(self, *treants, **kwargs):
        """Generate a Bundle from any number of Treants.

        :Arguments:
            *treants*
                list giving either Treants, Groups, or paths giving the
                directories of the state files for such objects in the
                filesystem

        :Keywords:
            *flatten* [NOT IMPLEMENTED]
                if ``True``, will recursively obtain members of any Groups;
                only Treants will be present in the bunch

        """
        self._backend = _BundleBackend()
        self._cache = dict()
        self._data = None

        self.add(*treants)

    def __repr__(self):
        return "<Bundle({})>".format(self._list())
