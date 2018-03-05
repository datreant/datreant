"""Trees and Leaves: filesystem manipulation interfaces for directories and
files.

"""
from __future__ import absolute_import

import os
from functools import reduce, total_ordering
from six import string_types

from collections import OrderedDict
import scandir
from pathlib2 import Path
from asciitree import LeftAligned

from .util import makedirs
from .manipulators import discover
from .rsync import rsync


@total_ordering
class Veg(object):

    def __init__(self, filepath):
        self._path = Path(os.path.abspath(filepath))

    def __fspath__(self):
        return str(self._path.absolute())

    def __str__(self):
        return str(self.path)

    def __eq__(self, other):
        try:
            return self.abspath == other.abspath
        except AttributeError:
            return NotImplemented

    def __lt__(self, other):
        try:
            return self.abspath < other.abspath
        except AttributeError:
            return NotImplemented

    def __hash__(self):
        return hash(self.abspath)

    @property
    def exists(self):
        """Check existence of this path in filesystem.

        """
        return self.path.exists()

    @property
    def path(self):
        """Filesystem path as a :class:`pathlib2.Path`.

        """
        return self._path

    @property
    def abspath(self):
        """Absolute path.

        """
        return str(self.path.absolute())

    @property
    def relpath(self):
        """Relative path from current working directory.

        """
        return os.path.relpath(str(self.path))

    @property
    def parent(self):
        """Parent directory for this path.

        """
        return Tree(str(self.path.parent))

    @property
    def name(self):
        """Basename for this path.

        """
        return os.path.basename(os.path.abspath(self.abspath))


class Leaf(Veg):
    """A file in the filesystem."""

    def __init__(self, filepath):
        if os.path.isdir(filepath):
            raise ValueError("'{}' is an existing directory; "
                             "a Leaf must be a file".format(filepath))

        self._path = Path(os.path.abspath(filepath))

    def __repr__(self):
        return "<Leaf: '{}'>".format(self.relpath)

    def makedirs(self):
        """Make all directories along path that do not currently exist.

        Returns
        -------
        leaf : Leaf
            this leaf

        """
        makedirs(os.path.dirname(str(self.path)), exist_ok=True)

        return self

    def touch(self):
        """Make file if it doesn't exist.

        """
        self.makedirs()
        self.path.touch()

    def make(self):
        """Make the file if it doesn't exit. Equivalent to :meth:`touch`.

        Returns
        -------
        leaf : Leaf
            this leaf

        """
        self.touch()
        return self

    def read(self, size=None):
        """Read file, or up to `size` in bytes.

        Parameters
        ----------
        size : int
            extent of the file to read, in bytes

        """
        with open(self.abspath, 'r') as f:
            if size:
                out = f.read(size)
            else:
                out = f.read()
        return out


class Tree(Veg):
    """A directory."""
    def __init__(self, dirpath):
        if isinstance(dirpath, Tree):
            dirpath = dirpath.abspath

        elif os.path.isfile(dirpath):
            raise ValueError("'{}' is an existing file; "
                             "a Tree must be a directory".format(dirpath))

        self._path = Path(os.path.abspath(dirpath))

    def __repr__(self):
        return "<Tree: '{}'>".format(self.relpath)

    def __contains__(self, item):
        """Returns True if given Tree, Leaf, or plain path resolves as being
        within this Tree.

        The given path need not exist in the filesystem.

        """
        if isinstance(item, (Tree, Leaf)):
            return str(self) in str(item)
        elif isinstance(item, string_types):
            return str(self) in os.path.abspath(item)
        else:
            raise TypeError("Item must be a Tree, Leaf, or plain path")

    def __getitem__(self, path):
        """Get trees or leaves in this tree.

        Parameters
        ----------
        path : str
            Path to desired tree or leaf, relative to this tree; a list of
            paths will give a View of the resulting Trees and Leaves

        Returns
        -------
        result : Tree or Leaf, View
            Tree or Leaf corresponding to given path; if multiple paths
            given in a list, a View of the resulting Trees or Leaves is
            returned

        """
        from .collections import View

        def filt(path):
            fullpath = os.path.abspath(os.path.join(self.abspath, path))

            if (os.path.isdir(fullpath) or path.endswith(os.sep) or
                    (fullpath in self.abspath)):
                return Tree(fullpath)
            else:
                return Leaf(fullpath)

        if isinstance(path, list):
            outview = []
            for item in path:
                outview.append(filt(item))

            return View(outview)
        elif isinstance(path, string_types):
            return filt(path)
        else:
            raise ValueError('Must use a path or a list of paths')

    @property
    def loc(self):
        """Get Tree/Leaf at relative `path`.

        Use with getitem syntax, e.g. ``.loc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If directory/file does not exist at the given path, then whether a Tree
        or Leaf is given is determined by the path semantics, i.e. a trailing
        separator ("/").

        Using e.g. ``Tree.loc['some name']`` is equivalent to doing
        ``Tree['some name']``. ``.loc`` is included for parity with ``View``
        and ``Bundle`` API semantics.

        """
        if not hasattr(self, "_loc"):
            self._loc = _Loc(self)

        return self._loc

    @property
    def treeloc(self):
        """Get Tree at relative `path`.

        Use with getitem syntax, e.g. ``.treeloc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If the given path resolves to an existing file, then a ``ValueError``
        will be raised.

        """
        if not hasattr(self, "_treeloc"):
            self._treeloc = _TreeLoc(self)

        return self._treeloc

    @property
    def leafloc(self):
        """Get Leaf at relative `path`.

        Use with getitem syntax, e.g. ``.treeloc['some name']``

        Allowed inputs are:
        - A single name
        - A list or array of names

        If the given path resolves to an existing directory, then a
        ``ValueError`` will be raised.

        """
        if not hasattr(self, "_leafloc"):
            self._leafloc = _LeafLoc(self)

        return self._leafloc

    @property
    def abspath(self):
        """Absolute path of ``self.path``."""
        return str(self.path.absolute()) + os.sep

    @property
    def relpath(self):
        """Relative path of ``self.path`` from current working directory."""
        return os.path.relpath(str(self.path)) + os.sep

    @property
    def parent(self):
        """Parent directory for this path."""
        return Tree(str(self.path.parent))

    def leaves(self, hidden=False):
        """Return a View of the files in this Tree.

        Parameters
        ----------
        hidden : bool
            If True, include hidden files.

        Returns
        -------
        View
            A View with files in this Tree as members.

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        for root, dirs, files in scandir.walk(self.abspath):
            if hidden:
                out = [Leaf(os.path.join(root, f)) for f in files]
            else:
                out = [Leaf(os.path.join(root, f)) for f in files
                       if f[0] != os.extsep]
            break

        out.sort()
        return View(out)

    def trees(self, hidden=False):
        """Return a View of the directories in this Tree.

        Parameters
        ----------
        hidden : bool
            If True, include hidden directories.

        Returns
        -------
        View
            A View with directories in this Tree as members.

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        for root, dirs, files in scandir.walk(self.abspath):
            if hidden:
                out = [Tree(os.path.join(root, d))
                       for d in dirs]
            else:
                out = [Tree(os.path.join(root, d))
                       for d in dirs if d[0] != os.extsep]
            break

        out.sort()
        return View(out)

    def children(self, hidden=False):
        """Return a View of all files and directories in this Tree.

        Parameters
        ----------
        hidden : bool
            If True, include hidden files and directories.

        Returns
        -------
        View
            A View with files and directories in this Tree as members.

        """
        from .collections import View
        return View((self.trees(hidden=hidden) +
                     self.leaves(hidden=hidden)))

    def glob(self, pattern):
        """Return a View of all child Leaves and Trees matching given globbing
        pattern.

        Parameters
        ----------
        pattern
            globbing pattern to match files and directories with

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        out = []
        for item in self.path.glob(pattern):
            out.append(self[str(item)])

        out.sort()

        return View(out)

    def walk(self, topdown=True, onerror=None, followlinks=False):
        """Walk through the contents of the tree.

        For each directory in the tree (including the root itself), yields a
        3-tuple (dirpath, dirnames, filenames).

        Parameters
        ----------
        topdown : Boolean, optional
            If False, walks directories from the bottom-up.
        onerror : function, optional
            Optional function to be called on error.
        followlinks : Boolean, optional
            If False, excludes symbolic file links.

        Returns
        -------
        generator
            Wrapped `scandir.walk()` generator yielding `datreant` objects

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        for root, dirs, files in os.walk(self.abspath, topdown=topdown,
                                         onerror=onerror,
                                         followlinks=followlinks):

            # wrap results in datreant objects
            r_tree = Tree(root)
            trees = r_tree[dirs]
            leaves = r_tree[files]

            yield r_tree, trees, leaves

    def draw(self, depth=None, hidden=False):
        """Print an ASCII-fied visual of the tree.

        Parameters
        ----------
        depth : int, optional
            Maximum directory depth to display. ``None`` indicates no limit.
        hidden : bool
            If True, show hidden files and directories.

        """
        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        tree = OrderedDict()
        rootdir = self.abspath.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in scandir.walk(rootdir):

            # sort files and directories so they show up in output sorted
            dirs.sort()
            files.sort()

            folders = ["{}/".format(x) for x in path[start:].split(os.sep)]
            parent = reduce(dict.get, folders[:-1], tree)

            # depth handling
            if depth and len(folders) == depth+1:
                parent[folders[-1]] = {}
                continue
            elif depth and len(folders) > depth+1:
                continue

            # filter out hidden files and directories, if desired
            if not hidden:
                outfiles = [file for file in files if file[0] != os.extsep]
                hidden_dirs = [d for d in dirs if d[0] == os.extsep]
                for d in hidden_dirs:
                    dirs.remove(d)
            else:
                outfiles = files

            subdir = OrderedDict.fromkeys(outfiles, {})
            parent[folders[-1]] = subdir

        tr = LeftAligned()
        print(tr(tree))

    def makedirs(self):
        """Make all directories along path that do not currently exist.

        Returns
        -------
        Tree
            This Tree.

        """
        makedirs(str(self.path), exist_ok=True)

        return self

    def make(self):
        """Make the directory if it doesn't exist. Equivalent to :meth:`makedirs`.

        Returns
        -------
        Tree
            This Tree.

        """
        self.makedirs()

        return self

    def sync(self, other, mode='upload', compress=True, checksum=True,
             backup=False, dry=False, include=None, exclude=None,
             overwrite=False, rsync_path='/usr/bin/rsync'):
        """Synchronize directories using rsync.

        Parameters
        ----------
        other: str or Tree
            Other end of the sync, can be either a path or another Tree.
        mode: str
            Either ``"upload"`` if uploading to  `other`, or ``"download"`` if
            downloading from `other`


        The other options are described in the
        :py:func:`datreant.rsync.rsync` documentation.

        """
        if isinstance(other, Tree):
            other = other.abspath

        if mode == 'download':
            source = other
            dest = self.abspath
        elif mode == 'upload':
            source = self.abspath
            dest = other
        else:
            raise ValueError("Sync mode can be only 'upload' or 'download'.")
        # Here we do some massaging before passing to the rsync function
        return rsync(source, dest, compress=compress, backup=backup,
                     dry=dry, include=include, checksum=checksum,
                     overwrite=overwrite, exclude=exclude,
                     rsync_path=rsync_path)


class _Loc(object):
    """Path accessor for Trees."""

    def __init__(self, tree):
        self._tree = tree

    def __getitem__(self, path):
        """Get Tree/Leaf at `path` relative to attached Tree.

        """
        return self._tree[path]


class _TreeLoc(_Loc):
    """Tree accessor for Trees."""

    def __getitem__(self, path):
        """Get Tree at `path` relative to attached Tree.

        """
        return Tree(self._tree[path].abspath)


class _LeafLoc(_Loc):
    """Leaf accessor for Trees."""

    def __getitem__(self, path):
        """Get Leaf at `path` relative to attached Tree.

        """
        return Leaf(self._tree[path].abspath)
