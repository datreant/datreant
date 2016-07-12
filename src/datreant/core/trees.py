"""Trees and Leaves: filesystem manipulation interfaces for directories and
files.

"""

import os
from functools import reduce, total_ordering
from six import string_types

import scandir
from pathlib2 import Path
from asciitree import LeftAligned

from .util import makedirs
from .manipulators import discover
from .rsync import rsync
from . import _TREELIMBS


@total_ordering
class Veg(object):

    _classlimbs = set()
    _limbs = set()

    def __init__(self, filepath):
        self._path = Path(os.path.abspath(filepath))

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
        if isinstance(self, Tree):
            limbs = self._classlimbs | self._limbs
        else:
            limbs = None

        return Tree(str(self.path.parent), limbs=limbs)

    @property
    def name(self):
        """Basename for this path.

        """
        return os.path.basename(os.path.abspath(self.abspath))

    @property
    def limbs(self):
        """A set giving the names of this object's attached limbs.

        """
        return self._classlimbs | self._limbs


class Leaf(Veg):
    """A file in the filesystem.

    """

    def __init__(self, filepath):
        if os.path.isdir(filepath):
            raise ValueError("'{}' is an existing directory; "
                             "a Leaf must be a file".format(filepath))

        self._path = Path(os.path.abspath(filepath))

    def __repr__(self):
        return "<Leaf: '{}'>".format(self.relpath)

    def makedirs(self):
        """Make all directories along path that do not currently exist.

        :Returns:
            *leaf*
                this leaf

        """
        makedirs(os.path.dirname(str(self.path)))

        return self

    def touch(self):
        """Make file if it doesn't exist.

        """
        self.makedirs()
        self.path.touch()

    def make(self):
        """Make the file if it doesn't exit. Equivalent to :meth:`touch`.

        :Returns:
            *leaf*
                this leaf

        """
        self.touch()
        return self

    def read(self, size=None):
        """Read file, or up to `size` in bytes.

        :Arguments:
            *size*
                extent of the file to read, in bytes

        """
        with open(self.abspath, 'r') as f:
            if size:
                out = f.read(size)
            else:
                out = f.read()
        return out


class Tree(Veg):
    """A directory.

    """
    def __init__(self, dirpath, limbs=None):
        if os.path.isfile(dirpath):
            raise ValueError("'{}' is an existing file; "
                             "a Tree must be a directory".format(dirpath))

        self._path = Path(os.path.abspath(dirpath))

        if limbs:
            self.attach(*limbs)

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
                limbs = self._classlimbs | self._limbs
                return Tree(fullpath, limbs=limbs)
            else:
                return Leaf(fullpath)

        if isinstance(path, list):
            outview = []
            for item in path:
                outview.append(filt(item))

            return View(outview, limbs=self.limbs)
        elif isinstance(path, string_types):
            return filt(path)
        else:
            raise ValueError('Must use a path or a list of paths')

    @classmethod
    def _attach_limb_class(cls, limb):
        """Attach a limb to the class.

        """
        # property definition
        def getter(self):
            if not hasattr(self, "_"+limb._name):
                setattr(self, "_"+limb._name, limb(self))
            return getattr(self, "_"+limb._name)

        try:
            setter = limb._setter
        except AttributeError:
            setter = None

        # set the property
        setattr(cls, limb._name,
                property(getter, setter, None, limb.__doc__))

        if limb._name in _TREELIMBS:
            cls._classlimbs.add(limb._name)

    def _attach_limb(self, limb):
        """Attach a limb.

        """
        try:
            setattr(self, limb._name, limb(self))
        except AttributeError:
            pass

        if limb._name in _TREELIMBS:
            self._limbs.add(limb._name)

    def attach(self, *limbname):
        """Attach limbs by name to this Tree.

        """
        for ln in limbname:
            try:
                limb = _TREELIMBS[ln]
            except KeyError:
                raise KeyError("No such limb '{}'".format(ln))
            else:
                self._attach_limb(limb)

    @property
    def loc(self):
        """Get Tree/Leaf at `path` relative to Tree.

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
    def abspath(self):
        """Absolute path of ``self.path``.

        """
        return str(self.path.absolute()) + os.sep

    @property
    def relpath(self):
        """Relative path of ``self.path`` from current working directory.

        """
        return os.path.relpath(str(self.path)) + os.sep

    @property
    def leaves(self):
        """A View of the files in this Tree.

        Hidden files are not included.

        """
        from .collections import View

        if self.exists:
            for root, dirs, files in scandir.walk(self.abspath):
                # remove hidden files
                out = [Leaf(os.path.join(root, f)) for f in files
                       if f[0] != os.extsep]
                break

            out.sort()
            return View(out, limbs=self.limbs)
        else:
            raise OSError("Tree doesn't exist in the filesystem")

    @property
    def trees(self):
        """A View of the directories in this Tree.

        Hidden directories are not included.

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        for root, dirs, files in scandir.walk(self.abspath):
            # remove hidden directories
            out = [Tree(os.path.join(root, d), limbs=self.limbs) for d in dirs
                   if d[0] != os.extsep]
            break

        out.sort()
        return View(out, limbs=self.limbs)

    @property
    def hidden(self):
        """A View of the hidden files and directories in this Tree.

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")
        for root, dirs, files in scandir.walk(self.abspath):
            outdirs = [Tree(os.path.join(root, d), limbs=self.limbs)
                       for d in dirs if d[0] == os.extsep]
            outdirs.sort()

            outfiles = [Leaf(os.path.join(root, f)) for f in files
                        if f[0] == os.extsep]
            outfiles.sort()

            # want directories then files
            out = outdirs + outfiles
            break

        return View(out, limbs=self.limbs)

    @property
    def children(self):
        """A View of all files and directories in this Tree.

        Includes hidden files and directories.

        """
        from .collections import View
        return View(self.trees + self.leaves + self.hidden, limbs=self.limbs)

    discover = discover

    @property
    def treants(self):
        """Bundle of all Treants found within this Tree.

        This does not return a Treant for a bare state file found within this
        Tree. In effect this gives the same result as ``Bundle(self.trees)``.

        """
        from .collections import Bundle
        return Bundle(self.trees + self.hidden.membertrees, limbs=self.limbs)

    def glob(self, pattern):
        """Return a View of all child Leaves and Trees matching given globbing
        pattern.

        :Arguments:
            *pattern*
               globbing pattern to match files and directories with

        """
        from .collections import View

        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        out = []
        for item in self.path.glob(pattern):
            out.append(self[str(item)])

        out.sort()

        return View(out, limbs=self.limbs)

    def draw(self, depth=None, hidden=False):
        """Print an ASCII-fied visual of the tree.

        Parameters
        ----------
        depth : int
            Maximum directory depth to display. ``None`` indicates no limit.
        hidden : bool
            If False, do not show hidden files; hidden directories are still
            shown if they contain non-hidden files or directories.

        """
        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")

        tree = {}
        rootdir = self.abspath.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in scandir.walk(rootdir):
            folders = ["{}/".format(x) for x in path[start:].split(os.sep)]

            parent = reduce(dict.get, folders[:-1], tree)

            if depth and len(folders) == depth+1:
                parent[folders[-1]] = {}
                continue
            elif depth and len(folders) > depth+1:
                continue

            # filter out hidden files, if desired
            if not hidden:
                outfiles = [file for file in files if file[0] != os.extsep]
            else:
                outfiles = files

            subdir = dict.fromkeys(outfiles, {})
            parent[folders[-1]] = subdir

        tr = LeftAligned()
        print(tr(tree))

    def makedirs(self):
        """Make all directories along path that do not currently exist.

        :Returns:
            *tree*
                this tree

        """
        makedirs(str(self.path))

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
        :py:func:`datreant.core.rsync.rsync` documentation.

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
    """Subtree accessor for Trees."""

    def __init__(self, tree):
        self._tree = tree

    def __getitem__(self, path):
        """Get Tree/Leaf at `path` relative to attached Tree.

        """
        return self._tree[path]
