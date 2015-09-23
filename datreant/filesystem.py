"""
Functions and classes for finding Treants in the filesystem.

"""
import os
import sys
import glob
import time

import scandir

from datreant import persistence
import datreant


def statefilename(treanttype, uuid, ext):
    """Return state file name given the type of treant and its uuid.

    """
    return "{}.{}{}".format(treanttype, uuid, ext)


def glob_treant(treant):
    """Given a Treant's directory, get its state file.

    Since state file names contain uuids, they vary. Multiple state files may
    therefore be present in a single directory. All will be returned as a list.

    :Arguments:
        *treant*
            directory containing a state file

    :Returns:
        *treantfile*
            list giving absolute paths of state files found
            in directory
    """
    fileglob = list()
    for treanttype in datreant._treants:
        for backend in datreant._treants[treanttype]._backends:
            extension = datreant._treants[treanttype]._backends[backend][0]
            fileglob.extend(
                glob.glob(os.path.join(
                    treant,
                    '{}.*{}'.format(treanttype, extension))))

    paths = [os.path.abspath(x) for x in fileglob]
    return paths


def path2treant(*paths):
    """Return Treants from directories or full paths containing Treant
        state files.

    *Note*: If there are multiple state files in a given directory, Treants
            will be returned for each.

    :Arguments:
        *paths*
            directories containing state files or full paths to state files to
            load Treants from; if ``None`` is an element, then ``None``
            returned in output list

    :Returns:
        *treants*
            list of Treants obtained from directories; ``None`` as an
            element indicates that ``None`` was present in the list of paths

    """
    treants = list()
    for path in paths:
        if path is None:
            treants.append(None)
        elif os.path.isdir(path):
            files = glob_treant(path)
            for item in files:
                basename = os.path.basename(item)
                for treanttype in datreant._treants:
                    if treanttype in basename:
                        treants.append(datreant._treants[treanttype](item))
        elif os.path.exists(path):
            basename = os.path.basename(path)
            for treanttype in datreant._treants:
                if treanttype in basename:
                    treants.append(datreant._treants[treanttype](path))

    return treants


class Foxhound(object):
    """Locator for Treants.

    This object is used by Treants to find Treants, even when they are no
    longer in their last known location.

    Groups and Bundles uses this class to find members that have moved. All
    TreantFiles use this class to find their file on disk when it moves.

    """
    def __init__(self, caller, uuids, basedirs, coordinators=None, timeout=10):
        """Generate a Foxhound to track down Treants.

        :Arguments:
            *caller*
                object that summoned the Foxhound; needed to make sense
                of some path types, as well as for automated context
                for conducting the search
            *uuids*
                list of unique identifiers of Treants to find
            *basedirs*
                dict of basedirs to start searching around; keys may be
                'abspath' or 'relCont', and values should be lists of paths

        :Keywords:
            *coordinators*
                list of Coordinators to consult; if ``None``, involve no
                Coordinators
            *timeout*
                maximum time, in seconds, the Foxhound will spend fetching.

        """
        self.caller = caller
        self.uuids = uuids
        self.basedirs = basedirs
        self.coordinators = coordinators

        self.timeout = timeout

        # once found: uuids as keys, absolute paths as values
        self.treants = dict()

    def fetch(self, as_treants=True):
        """Find the Treants.

        :Keywords:
            *as_treants*
                if ``True``, return Treant instances instead of absolute
                paths to state files

        :Returns:
            *results*
                dictionary giving Treant uuids as keys and absolute paths to
                their state files as values; ``None`` as a value indicates
                that no state file could be found. Returns Treant instances
                instead of paths for *as_treants* == True.

        """
        from datreant.aggregators import Members
        from datreant.collections import Bundle

        if isinstance(self.caller, Members):
            results = self._find_Group_members()
        elif isinstance(self.caller, Bundle):
            results = self._find_Bundle_members()

        if as_treants:
            conts = path2treant(*results.values())
            results = {x: y for x, y in zip(results.keys(), conts)}

        return results

    def _check_basedirs(self):
        """Check last-known locations for Treants.

        :Returns:
            *results*
                dictionary giving Treant uuids as keys and absolute paths to
                their state files as values; ``None`` as a value indicates
                that no state file could be found.
        """
        # initialize output dictionary with None
        outpaths = {x: y for x, y in zip(self.uuids, [None]*len(self.uuids))}

        uuids = [x for x in outpaths if not outpaths[x]]
        if 'abspath' in self.basedirs:
            for path in self.basedirs['abspath']:
                found = []
                for uuid in uuids:
                    candidate = glob.glob(
                            os.path.join(path, '*.{}.*'.format(uuid)))

                    if candidate:
                        outpaths[uuid] = os.path.abspath(candidate[0])
                        found.append(uuid)

                for item in found:
                    uuids.remove(item)

        if 'relCont' in self.basedirs:
            # get uuids for which paths haven't been found
            for path in self.basedirs['relCont']:
                found = []
                for uuid in uuids:
                    candidate = glob.glob(
                            os.path.join(
                                self.caller._backend.get_location(),
                                path, '*.{}.*'.format(uuid)))

                    if candidate:
                        outpaths[uuid] = os.path.abspath(candidate[0])
                        found.append(uuid)

                for item in found:
                    uuids.remove(item)

        return outpaths

    def _downward_search(self, path):
        """Check for Treants downward from specified path.

        :Arguments:
            *path*
                path to begin downward search from

        """
        pass

    def _outward_search(self, path):
        pass

    def _consult_Coordinators(self):
        pass

    def _find_TreantFile(self):
        """Find Treant for a TreantFile.

        If a Treant's state file is moved by another process while a
        Treant instance exists, then the TreantFile instance must
        find its state file when it discovers it has gone missing. The
        Foxhound begins by searching downward from the Treant's previous
        location, with subsequent downward searches proceeding from the parent
        directory. This process continues until either the state file is found,
        the filesystem is exhaustively searched, or the Foxhound times out.

        """
        pass

    def _find_Group_members(self):
        """Find Treants that are members of a Group.

        For finding Group members, the Foxhound begins by looking for
        Treants among the paths it was given. Treants that can't be found
        are then searched for starting downward from the Group's location, with
        subsequent downward searches proceeding from the parent directory.
        This process continues until either all members are found, the
        filesystem is exhaustively searched, or the Foxhound times out.

        :Returns:
            *outpaths*
                dictionary giving Treant uuids as keys and absolute paths to
                their state files as values; ``None`` as a value indicates
                that no state file could be found.

        """
        # search last-known locations
        outpaths = self._check_basedirs()

        # get current time
        currtime = time.time()

        # walk downwards on an upward path through filesystem from the Group's
        # basedir
        uuids = [x for x in outpaths if not outpaths[x]]
        path = self.caller._backend.get_location()
        prev = None
        timedout = False
        while prev != path and uuids and not timedout:

            top = True
            for root, dirs, files in scandir.walk(path):
                # if search runs over timeout, call it off
                if ((time.time() - currtime) > self.timeout):
                    self.caller._logger.info(
                            "Search for missing members timed" +
                            " out at {}".format(self.timeout) +
                            " seconds.")
                    timedout = True
                    break

                # if we've found everything, finish
                if not uuids:
                    break

                found = []
                # no need to visit already-visited tree
                if top and prev:
                    dirs.remove(os.path.basename(prev))
                    top = False

                for uuid in uuids:
                    candidate = [os.path.join(root, x)
                                 for x in files if (uuid in x and x[0] != '.')]

                    if candidate:
                        outpaths[uuid] = os.path.abspath(candidate[0])
                        found.append(uuid)

                for item in found:
                    uuids.remove(item)

            prev = path
            path = os.path.split(path)[0]

        # TODO: post-check? Since Groups know the treanttypes of their
        # members, should we compare these to what is in outpaths?

        return outpaths

    def _find_Bundle_members(self):
        """Find Treants that are members of a Bundle.

        For finding Bundle members, the Foxhound begins by looking for
        Treants among the paths it was given. Treants that can't be found
        are then searched for starting downward from the current working
        directory with subsequent downward searches proceeding from the parent
        directory. This process continues until either all members are found,
        the filesystem is exhaustively searched, or the Foxhound times out.

        :Returns:
            *outpaths*
                dictionary giving Treant uuids as keys and absolute paths to
                their state files as values; ``None`` as a value indicates
                that no state file could be found.

        """
        # search last-known locations
        outpaths = self._check_basedirs()

        # get current time
        currtime = time.time()

        # walk downwards on an upward trajectory through filesystem from the
        # current working directory
        uuids = [x for x in outpaths if not outpaths[x]]
        path = os.path.abspath(os.curdir)
        prev = None
        while prev != path and uuids:

            # if search runs over timeout, call it off
            if ((time.time() - currtime) > self.timeout):
                self.caller._logger.info("Search for missing members timed" +
                                         " out at {}".format(self.timeout) +
                                         " seconds.")
                break

            top = True
            for root, dirs, files in scandir.walk(path):
                found = []
                # no need to visit already-visited tree
                if top and prev:
                    dirs.remove(os.path.basename(prev))
                    top = False

                for uuid in uuids:
                    candidate = [os.path.join(root, x)
                                 for x in files if (uuid in x and x[0] != '.')]

                    if candidate:
                        outpaths[uuid] = os.path.abspath(candidate[0])
                        found.append(uuid)

                for item in found:
                    uuids.remove(item)

            prev = path
            path = os.path.split(path)[0]

        # TODO: post-check? Since Bundles know the treanttypes of their
        # members, should we compare these to what is in outpaths?

        return outpaths

    def _find_Coordinator_members(self):
        pass

    def discover(self, path):
        pass

    # OLD
    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Treant looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The directory
        containing the first such file found will be returned. None is returned
        if no such files found.

        :Keywords:
            *startdir*
                directory from which to begin upward search; default is
                Treant basedir

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
