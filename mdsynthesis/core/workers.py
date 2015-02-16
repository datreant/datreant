"""
Under-the-hood classes. Mostly a grab-bag of needed functionality,
and possibly a bit messy.

"""
import os, sys
import glob

import aggregators
import persistence
import mdsynthesis as mds

def statefilename(containertype, uuid):
    """Return state file name given the type of container and its uuid.

    """
    return "{}.{}.{}".format(containertype, uuid, persistence.statefile_ext)

def glob_containerfile(container):
    """Given a Container's directory, get its state file.

    Since state file names contain uuids, they vary.

    :Arguments:
        *container*
            directory containing a state file

    :Returns:
        *containerfile*
            list giving absolute paths of state files found
            in directory
    """


class Utilities(object):
    """Mixin with a few commonly-used standalone methods.

    If it needs to be used in more than one unrelated class, it goes here.

    """
    def _path2container(self, *directories):
        """Return Containers from directories containing Container state files.

        :Arguments:
            *directories*
                directories containing state files to be loaded from
    
        :Returns:
            list of Containers obtained from directories; returns ``None`` for
            paths that didn't yield a Container

        """
        containers = []
        for directory in directories:
            if os.path.exists(os.path.join(directory, persistence.simfile)):
                containers.append(mds.Sim(directory))
            elif os.path.exists(os.path.join(directory, persistence.groupfile)):
                containers.append(mds.Group(directory))
            else:
                containers.append(None)
    
        return containers

class Foxhound(object):
    """Locator for Containers that have gone missing.

    This object is used by Containers to find Containers when they are no longer
    in their last known location.

    """
    def __init__(self, uuids, locations, coordinators=None):
        """Generate a foxhound to track down Containers.

        :Arguments:
            *uuids*
                list of unique identifiers of Containers to find
            *locations*
                list of locations to start searching around; these need not
                be in any particular order, nor must their number match that
                of the given uuids
            *coordinators*
                list of Coordinators to consult; if ``None``, involve no Coordinators

        """
        self.uuids = uuids
        self.locations = locations
        self.coordinators = coordinators

        # once found: uuids as keys, absolute paths as values
        self.containers = dict()

    def _downward_search(self, path):
        """Check for Containers downward from specified path.

        :Arguments:
            *path*
                path to begin downward search from

        """
        pass

    def _outward_search(self, path):
        pass

    def _consult_Coordinators(self):
        pass

    def find_Group_members(self):
        pass

    def find_Coordinator_members(self):
        pass

    def discover(self, path):
        pass

    #OLD
    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The directory containing
        the first such file found will be returned. None is returned if no such
        files found.

        :Keywords:
            *startdir*
                directory from which to begin upward search; default is
                Container basedir

        :Returns:
            *database*
                directory of located Database; if no Database found, is None
        
        """
        startdir = kwargs.pop('startdir', None)
        
        if not startdir:
            startdir = self.metadata['basedir']

        # search upward for a database
        startdir = os.path.abspath(startdir)
        directory = startdir
        found = False
        
        self._logger.info("Beginning search for database from {}".format(directory))

        while (directory != '/') and (not found):
            directory, tail = os.path.split(directory)
            candidates = glob.glob(os.path.join(directory, self._databasefile))
            
            if candidates:
                self._logger.info("Database candidate located: {}".format(candidates[0]))
                basedir = os.path.dirname(candidates[0])
                db = Database.Database(basedir)
                found = db._handshake()
        
        if not found:
            self._logger.warning("No database found!")
            basedir = None

        return basedir

class Bundle(Utilities):
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
                cont = self._path2container(container)[0]
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
