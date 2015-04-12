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
            self._backend.add_member(container.uuid, container.name, container.containertype, container.location)

    def remove(self, *indices, **kwargs): 
        """Remove any number of members from the Group.
    
        :Arguments:
            *indices*
                the indices of the members to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        uuids = self._containerfile.get_members_uuid()
        if not kwargs.pop('all', False):
            uuids = [ uuids[x] for x in indices ]

        self._backend.del_member(*uuids)

#TODO: use an in-memory SQLite db instead of a list for this?
class _BundleBackend():
    """Backend class for Bundle. 
    
    Has same interface as Group-specific components of
    :class:`persistence.GroupFile`. Behaves practically like an in-memory
    version of a state-file, but with only the components needed for the
    Bundle.

    """

    def __init__(self):
        # our table will be a list of dicts
        self.table = list()

    @File._write_state
    def add_member(self, uuid, name, containertype, location):
        """Add a member to the Group.

        If the member is already present, its location will be updated with
        the given location.

        :Arguments:
            *uuid*
                the uuid of the new member
            *containertype*
                the container type of the new member (Sim or Group)
            *location*
                location of the new member in the filesystem
    
        """
        # check if uuid already present
        index = [ self.table.index(item) for item in self.table if item['uuid'] == uuid ]
        if index:
            # if present, update location
            self.table[index]['abspath'] = os.path.abspath(location)
        else:
            newmem = {'uuid': uuid,
                      'name': name,
                      'containertype': containertype,
                      'abspath': os.path.abspath(location)}
                      
            self.table.append(newmem)

    @File._write_state
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
            self.table = list()
        else:
            # remove redundant uuids from given list if present
            uuids = set([ str(uid) for uid in uuid ])

            # remove matching elements
            matches = list()
            for element in self.table:
                for uuid in uuids:
                    if (element['uuid'] == uuid):
                        matches.append(element)

            for element in matches:
                self.table.remove(element)

    @File._read_state
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
        # check if uuid present
        index = [ self.table.index(item) for item in self.table if item['uuid'] == uuid ]
        if index:
            memberinfo = self.table[index]
        else:
            self.logger.info('No such member.')
            memberinfo = None

        return memberinfo

    @File._read_state
    def get_members_uuid(self):
        """List uuid for each member.

        :Returns:
            *uuids*
                list giving uuids of all members, in order

        """
        return [ x['uuid'] for x in self.table ]

    @File._read_state
    def get_members_name(self):
        """List name for each member.

        :Returns:
            *names*
                list giving names of all members, in order

        """
        return [ x['name'] for x in self.table ]

    @File._read_state
    def get_members_containertype(self):
        """List containertype for each member.

        :Returns:
            *containertypes*
                list giving containertypes of all members, in order

        """
        return [ x['containertype'] for x in self.table ]

    @File._read_state
    def get_members_location(self, path='abspath'):
        """List stored location for each member. 

        :Arguments:
            *path*
                type of paths to return; only absolute paths (abspath)
                are stored for Bundles

        :Returns:
            *locations*
                list giving locations of all members, in order

        """
        return [ x[path] for x in table.iterrows() ]

class Bundle(object):
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

        self.add(*containers)

    #TODO: make more efficient by using ordered dict for storage?
    def add(self, *containers):
        outconts = list()
        for container in containers:
            if isinstance(container, list):
                self.add(*container)
            elif isinstance(container, mds.containers.Container):
                uuid = container.uuid
                if not (uuid in self._uuids):
                    outconts.append(container)
                    self._uuids.append(uuid)
            elif os.path.exists(container):
                conts = filesystem.path2container(container)
                for cont in conts:
                    uuid = cont.uuid
                    if not (uuid in self._uuids):
                        outconts.append(cont)
                        self._uuids.append(uuid)

        self._containers.extend(outconts)
    
    def _list(self):
        """Return list representation.
    
        """
        return list(self._containers)

    def __repr__(self):
        return "<Bundle({})>".format(self.list())

