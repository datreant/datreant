"""
The Bundle object is the primary manipulator for Containers in aggregate.
They are returned as queries to Groups, Coordinators, and other Bundles. They
offer convenience methods for dealing with many Containers at once.

"""
import os

import aggregators
import persistence
import filesystem
import mdsynthesis as mds

class _CollectionBase(object):
    """Common interface elements for ordered sets of Containers.

    :class:`aggregators.Members` and :class:`Bundle` both use this interface.

    """

    def __getitem__(self, index):
        """Get member corresponding to the given index or slice.
        
        """
        allrecords = self._backend.get_members()
        records = allrecords[index]
        uuids = records['uuid']

        if isinstance(uuids, basestring):
            member = None
            # if member already cached, use cached member
            if uuids in self._cache:
                member = self._cache[uuids]
            else:
                #TODO: STOPPED HERE

        elif isinstance(uuids, np.ndarray):
            member = list()
            for uuid, record in zip(uuids, records):
                # if member already cached, use cached member
                if uuid in self._cache:
                    member.append(self._cache[uuid])
                else:
                    newmember = None
                    for pathtype in self._backend.memberpaths:
                        # use full path to state file in case there are multiples, and to avoid
                        # loading a replacement (checks uuid)
                        path = os.path.join(record[pathtype], 
                                filesystem.statefilename(record['containertype'], record['uuid']))

                    if not newmember:
                        ind = list(allrecords['uuid']).index(uuid)
                        raise IOError("Could not find member {} (uuid: {}); re-add or remove it.".format(ind, uuid))

        raise IOError("Could not find member {} (uuid: {}); re-add or remove it.".format(index, uuids))

        return member

    def add(self, *containers):
        """Add any number of members to this collection.

        :Arguments:
            *containers*
                Sims and/or Groups to be added; may be a list of Sims and/or
                Groups; Sims or Groups can be given as either objects or paths
                to directories that contain object statefiles
        """
        outconts = list()
        for container in containers:
            if isinstance(container, (list, tuple, Bundle, aggregators.Members)):
                self.add(*container)
            elif isinstance(container, mds.Container):
                outconts.append(container)
            elif os.path.exists(container):
                cont = filesystem.path2container(container)
                for c in cont:
                    outconts.append(c)

        for container in outconts:
            self._backend.add_member(container.uuid, container.containertype, container.basedir)

    def remove(self, *indices, **kwargs): 
        """Remove any number of members from the Group.
    
        :Arguments:
            *indices*
                the indices of the members to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        uuids = self._backend.get_members_uuid()
        if not kwargs.pop('all', False):
            uuids = [ uuids[x] for x in indices ]

        self._backend.del_member(*uuids)

        # remove from cache
        for uuid in uuids:
            self._cache.pop(uuid, None)

    def containertypes(self):
        """Return a list of member containertypes.

        """
        return self._backend.get_members_containertype()

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

    def _list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        Missing members will be present in the list as ``None``. This method is not intended
        for user-level use.

        """
        members = self._backend.get_members()
        uuids = members['uuid']

        findlist = list()
        memberlist = list()

        for uuid in uuids:
            if uuid in self._cache:
                memberlist.append(self._cache[uuid])
            else:
                memberlist.append(None)
                findlist.append(uuid)

        # track down our non-cached containers
        foxhound = filesystem.Foxhound(self, findlist, locations, coordinators=self.coordinators)
        foundconts = foxhound.fetch()

        # insert found containers into output list
        for cont in foundconts:
            memberlist[memberlist.index(None)] = cont

        # add newly found containers to cache
        for uuid, cont in zip(findlist, foundconts):
            if cont:
                self._cache[uuid] = cont

        return memberlist

class _BundleBackend():
    """Backend class for Bundle. 
    
    Has same interface as Group-specific components of
    :class:`persistence.GroupFile`. Behaves practically like an in-memory
    version of a state-file, but with only the components needed for the
    Bundle.

    """

    def __init__(self):
        # our table will be a structured array matching the schema of the
        # GroupFile _Members Table
        self.table = None

    def _member2record(uuid, containertype, basedir):
        """Return a record array from a member's information.

        This method defines the scheme for the Bundle's record array.

        """
        return np.array((uuid, containertype, os.path.abspath(basedir)),
                        dtype={'names': ['uuid', 'containertype', 'abspath'],
                               'formats':['a{}'.format(persistence.uuidlength),
                                          'a{}'.format(persistence.namelength),
                                          'a{}'.format(persistence.pathlength)]})

    def add_member(self, uuid, containertype, basedir):
        """Add a member to the Group.

        If the member is already present, its location will be updated with
        the given location.

        :Arguments:
            *uuid*
                the uuid of the new member
            *containertype*
                the container type of the new member (Sim or Group)
            *basedir*
                basedir of the new member in the filesystem
    
        """
        if self.table is None:
            self.table = self._member2record(uuid, containertype, basedir)
        else:
            # check if uuid already present
            index = np.where(self.table['uuid'] == uuid)[0]
            if index:
                # if present, update location
                self.table[index[0]]['abspath'] = os.path.abspath(basedir)
            else:
                newmem = self._member2record(uuid, containertype, basedir)
                self.table = np.vstack(self.table, newmem)

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
            self.table = None
        else:
            # remove redundant uuids from given list if present
            uuids = set([ str(uid) for uid in uuid ])

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
            memberinfo = { x: memberinfo[x] for x in memberinfo.dtype.names }
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
                array giving containertype of each member, in order
        """
        return self.table['uuid']

    def get_members_containertype(self):
        """List containertype for each member.

        :Returns:
            *containertypes*
                array giving containertype of each member, in order
        """
        return self.table['containertype']

    def get_members_basedir(self):
        """List basedir for each member. 

        :Returns:
            *basedirs*
                structured array containing all paths to member basedirs
        """
        return self.table['abspath']

class Bundle(_CollectionBase):
    """Non-persistent Container for Sims and Groups.
    
    A Bundle is basically an indexable set. It is often used to return the
    results of a query on a Coordinator or a Group, but can be used on its
    own as well.

    """
    def __init__(self, *containers, **kwargs):
        """Generate a Bundle from any number of Containers.
    
        :Arguments:
            *containers*
                list giving either Sims, Groups, or paths giving the
                directories of the state files for such objects in the
                filesystem
    
        :Keywords:
            *flatten* [NOT IMPLEMENTED]
                if ``True``, will recursively obtain members of any Groups;
                only Sims will be present in the bunch 
         
        """
        self._backend = _BundleBackend()
        self._cache = dict()

        self.add(*containers)

    def __repr__(self):
        return "<Bundle({})>".format(self.list())

