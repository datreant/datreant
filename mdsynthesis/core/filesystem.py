"""
Functions and classes for finding Containers in the filesystem.

"""
import os
import sys
import glob

import aggregators
import persistence
import mdsynthesis as mds


def statefilename(containertype, uuid):
    """Return state file name given the type of container and its uuid.

    """
    return "{}.{}.{}".format(containertype, uuid, persistence.statefile_ext)


def glob_container(container):
    """Given a Container's directory, get its state file.

    Since state file names contain uuids, they vary. Multiple state files may
    therefore be present in a single directory. All will be returned as a list.

    :Arguments:
        *container*
            directory containing a state file

    :Returns:
        *containerfile*
            list giving absolute paths of state files found
            in directory
    """
    fileglob = list()
    for containertype in ('Container', 'Sim', 'Group'):
        fileglob.extend(
            glob.glob(os.path.join(container,
                                   '{}.*.h5'.format(containertype))))

    paths = [os.path.abspath(x) for x in fileglob]
    return paths


def path2container(*paths):
    """Return Containers from directories or full paths containing Container
        state files.

    *Note*: If there are multiple state files in a given directory, Containers
            will be returned for each.

    :Arguments:
        *paths*
            directories containing state files or full paths to state files to
            load Containers from

    :Returns:
        *containers*
            list of Containers obtained from directories

    """
    containers = list()
    for path in paths:
        if os.path.isdir(path):
            files = glob_container(path)
            for item in files:
                basename = os.path.basename(item)
                if 'Container' in basename:
                    containers.append(mds.Container(item))
                elif 'Sim' in basename:
                    containers.append(mds.Sim(item))
                elif 'Group' in basename:
                    containers.append(mds.Group(item))
        elif os.path.exists(path):
            basename = os.path.basename(path)
            if 'Container' in basename:
                containers.append(mds.Container(path))
            elif 'Sim' in basename:
                containers.append(mds.Sim(path))
            elif 'Group' in basename:
                containers.append(mds.Group(path))

    return containers


class Foxhound(object):
    """Locator for Containers that have gone missing.

    This object is used by Containers to find Containers when they are no
    longer in their last known location.

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
                list of Coordinators to consult; if ``None``, involve no
                Coordinators

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

    # OLD
    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The directory
        containing the first such file found will be returned. None is returned
        if no such files found.

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

        self._logger.info(
            "Beginning search for database from {}".format(directory))

        while (directory != '/') and (not found):
            directory, tail = os.path.split(directory)
            candidates = glob.glob(os.path.join(directory, self._databasefile))

            if candidates:
                self._logger.info(
                    "Database candidate located: {}".format(candidates[0]))
                basedir = os.path.dirname(candidates[0])
                db = Database.Database(basedir)
                found = db._handshake()

        if not found:
            self._logger.warning("No database found!")
            basedir = None

        return basedir
