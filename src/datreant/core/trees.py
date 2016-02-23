import os
from functools import reduce
from six import string_types

import scandir
from pathlib import Path
from asciitree import LeftAligned

from .util import makedirs
from . import _TREELIMBS


class BrushMixin(object):
    def __str__(self):
        return str(self.path)

    @property
    def exists(self):
        """Check existence of ``self.path`` in filesystem.

        """
        return self.path.exists()

    @property
    def path(self):
        """Filesystem as a :class:`pathlib.Path`.

        """
        return self._path

    @property
    def abspath(self):
        """Absolute path of ``self.path``.

        """
        return str(self.path.absolute())

    @property
    def relpath(self):
        """Relative path of ``self.path`` from current working directory.

        """
        return os.path.relpath(str(self.path))

    @property
    def parent(self):
        return Tree(str(self.path.parent))


class Leaf(BrushMixin):
    """A file in the filesystem.

    """

    def __init__(self, filepath):
        makedirs(os.path.dirname(filepath))
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

        """
        self.touch()


class Tree(BrushMixin):
    """A directory.

    """
    _classtreelimbs = []
    _treelimbs = []

    def __init__(self, dirpath, limbs=None):
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
            Path to desired tree or leaf, relative to this tree

        Returns
        -------
        Tree or Leaf

        """
        fullpath = os.path.join(self.abspath, path)

        if os.path.isdir(fullpath) or fullpath.endswith(os.sep):
            limbs = self._classtreelimbs + self._treelimbs
            return Tree(fullpath, limbs=limbs)
        else:
            return Leaf(fullpath)

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
            cls._classtreelimbs.append(limb._name)

    def _attach_limb(self, limb):
        """Attach a limb.

        """
        try:
            setattr(self, limb._name, limb(self))
        except AttributeError:
            pass

        if limb._name in _TREELIMBS:
            self._treelimbs.append(limb._name)

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
    def leaves(self):
        """A list of the file names in the directory.

        Hidden files are not included.

        """
        if self.exists:
            for root, dirs, files in scandir.walk(self.abspath):
                # remove hidden files
                out = [f for f in files if f[0] != os.extsep]
                break

            out.sort()
            return out
        else:
            raise OSError("Tree doesn't exist in the filesystem")

    @property
    def trees(self):
        """A list of the directories in the directory.

        Hidden directories are not included.

        """
        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")
        for root, dirs, files in scandir.walk(self.abspath):
            # remove hidden directories
            out = [d for d in dirs if d[0] != os.extsep]
            break

        out.sort()
        return out

    @property
    def hidden(self):
        """A list of the hidden files and directories in the directory.

        """
        if not self.exists:
            raise OSError("Tree doesn't exist in the filesystem")
        for root, dirs, files in scandir.walk(self.abspath):
            outdirs = [d for d in dirs if d[0] == os.extsep]
            outdirs.sort()

            outfiles = [f for f in files if f[0] == os.extsep]
            outfiles.sort()

            # want directories then files
            out = outdirs + outfiles
            break

        return out

    def draw(self, depth=None, hidden=False):
        """Print an asciified visual of the tree.

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

    @property
    def view(self):
        """Return contents of tree as a view.

        """
        pass

    def make(self):
        """Make the directory if it doesn't exit. Equivalent to :meth:`makedirs`.

        """
        self.makedirs()
