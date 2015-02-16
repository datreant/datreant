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
        self._containers = list()
        self._uuids = list()

        self.add(*containers)

    def add(self, *containers):
        outconts = list()
        for container in containers:
            if isinstance(container, list):
                self.add(*container)
            elif isinstance(container, mds.Sim) or isinstance(container, mds.Group):
                uuid = container._uuid
                if not (uuid in self._uuids):
                    outconts.append(container)
                    self._uuids.append(uuid)
            elif os.path.isdir(container):
                cont = filesystem.path2container(container)[0]
                if cont:
                    uuid = cont._uuid
                    if not (uuid in self._uuids):
                        outconts.append(cont)
                        self._uuids.append(uuid)

        self._containers.extend(outconts)
    
    def list(self):
        """Return list representation.
    
        """
        return list(self._containers)

    def __repr__(self):
        return "<Bundle({})>".format(self.list())
