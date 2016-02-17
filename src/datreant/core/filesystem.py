"""
Functions and classes for finding Treants in the filesystem.

"""
import os
import sys
import glob
import time

import scandir

from . import backends


def statefilename(treanttype, uuid):
    """Return state file name given the type of treant and its uuid.

    """
    return "{}.{}.{}".format(treanttype, uuid, 'json')


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
    from . import _TREANTS

    fileglob = list()
    for treanttype in _TREANTS:
        fileglob.extend(
            glob.glob(os.path.join(
                treant,
                '{}.*.{}'.format(treanttype, 'json'))))

    paths = [os.path.abspath(x) for x in fileglob]
    return paths


class TreantRec(object):
    def __init__(self, filepath):
        self.abspath = os.path.abspath(os.path.dirname(filepath))
        _basesplit = (os.path.basename(filepath)).split(os.extsep)

        self.treanttype = _basesplit[0]
        self.uuid = _basesplit[1]


def path2record(*paths):
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
    treantrecs = list()
    for path in paths:
        if path is None:
            treantrecs.append(None)
        elif os.path.isdir(path):
            files = glob_treant(path)
            for item in files:
                treantrecs.append(TreantRec(item))
        elif os.path.exists(path):
            treantrecs.append(TreantRec(path))

    return treantrecs


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
    from . import _TREANTS
    treants = list()
    for path in paths:
        if path is None:
            treants.append(None)
        elif os.path.isdir(path):
            files = glob_treant(path)
            for item in files:
                treanttype = os.path.basename(item).split(os.extsep)[0]
                treants.append(_TREANTS[treanttype](item))
        elif os.path.exists(path):
            treanttype = os.path.basename(path).split(os.extsep)[0]
            treants.append(_TREANTS[treanttype](path))

    return treants


class Foxhound(object):
    """Locator for Treants.

    This object is used by Treants to find Treants, even when they are no
    longer in their last known location.

    Groups and Bundles uses this class to find members that have moved. All
    TreantFiles use this class to find their file on disk when it moves.

    """
    def __init__(self, caller, uuids, paths, timeout=10):
        """Generate a Foxhound to track down Treants.

        :Arguments:
            *caller*
                object that summoned the Foxhound; needed to make sense
                of some path types, as well as for automated context
                for conducting the search
            *uuids*
                list of unique identifiers of Treants to find
            *paths*
                dict of directory paths to start searching around; keys may be
                'abspath' or 'relpath', and values should be lists of paths

        :Keywords:
            *timeout*
                maximum time, in seconds, the Foxhound will spend fetching.

        """
        self.caller = caller
        self.uuids = uuids
        self.paths = paths

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
        from .limbs import MemberBundle
        from .collections import Bundle

        if isinstance(self.caller, MemberBundle):
            results = self._find_Group_members()
        elif isinstance(self.caller, Bundle):
            results = self._find_Bundle_members()

        if as_treants:
            conts = path2treant(*results.values())
            results = {x: y for x, y in zip(results.keys(), conts)}

        return results

    def _check_paths(self):
        """Check last-known locations for Treants.

        :Returns:
            *results*
                dictionary giving Treant uuids as keys and absolute paths to
                their state files as values; ``None`` as a value indicates
                that no state file could be found.
        """
        # initialize output dictionary with None
        outpaths = dict.fromkeys(self.uuids)
        uuids = list(self.uuids)

        if 'abspath' in self.paths:
            for path in self.paths['abspath']:
                found = []
                candidates = glob.glob(
                        os.path.join(path, '*.*.json'))

                for candidate in candidates:
                    for uuid in uuids:
                        if uuid in candidate:
                            outpaths[uuid] = os.path.abspath(candidate)
                            found.append(uuid)
                            break

                for item in found:
                    uuids.remove(item)

        if 'relpath' in self.paths and uuids:
            import pdb
            pdb.set_trace()
            # get uuids for which paths haven't been found
            for path in self.paths['relpath']:
                found = []
                candidates = glob.glob(
                        os.path.join(path, '*.*.json'))

                for candidate in candidates:
                    for uuid in uuids:
                        if uuid in candidate:
                            outpaths[uuid] = os.path.abspath(candidate)
                            found.append(uuid)
                            break

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
        outpaths = self._check_paths()

        # get current time
        currtime = time.time()

        # walk downwards on an upward path through filesystem from the Group's
        # basedir
        uuids = [str(x) for x in outpaths if not outpaths[x]]
        path = self.caller._treant.location
        prev = None
        timedout = False
        while prev != path and uuids and not timedout:

            top = True
            for root, dirs, files in scandir.walk(path):
                # if search runs over timeout, call it off
                if ((time.time() - currtime) > self.timeout and
                        self.timeout is not None):

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
                    candidate = [os.path.join(root, x) for x in files
                                 if ((uuid in x) and
                                     ('.json' in x) and
                                     (x[0] != '.'))]

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
        outpaths = self._check_paths()

        # get current time
        currtime = time.time()

        # walk downwards on an upward trajectory through filesystem from the
        # current working directory
        uuids = [str(x) for x in outpaths if not outpaths[x]]
        path = os.path.abspath(os.curdir)
        prev = None
        while prev != path and uuids:
            import pdb
            pdb.set_trace()

            # if search runs over timeout, call it off
            if ((time.time() - currtime) > self.timeout and
                    self.timeout is not None):
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

    def discover(self, path):
        pass
