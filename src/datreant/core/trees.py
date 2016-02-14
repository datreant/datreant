import os
from six import string_types
import scandir
from pathlib import Path

from asciitree import LeftAligned

from .util import makedirs


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
        return str(self.path.relative_to(os.path.abspath('.')))


class Tree(BrushMixin):
    """A directory.

    """
    def __init__(self, dirpath):
        makedirs(dirpath)
        self._path = Path(os.path.abspath(dirpath))

    def __repr__(self):
        return "<Tree: '{}'>".format(self.relpath)

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

        if os.path.isdir(fullpath) or fullpath.endswith('/'):
            return Tree(fullpath)
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

        # set the property
        setattr(cls, limb._name,
                property(getter, None, None, limb.__doc__))

    def _attach_limb(self, limb):
        """Attach a limb.

        """
        try:
            setattr(self, limb._name, limb(self))
        except AttributeError:
            pass

    @property
    def leaves(self):
        """A list of the file names in the directory.

        Hidden files are not included.

        """
        for root, dirs, files in scandir.walk(self.abspath):
            # remove hidden files
            out = [f for f in files if f[0] != os.extsep]
            break

        out.sort()
        return out

    @property
    def trees(self):
        """A list of the directories in the directory.

        Hidden directories are not included.

        """
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
        print tr(tree)


class Leaf(BrushMixin):
    """A file in the filesystem.

    """

    def __init__(self, filepath):
        makedirs(os.path.dirname(filepath))
        self._path = Path(os.path.abspath(filepath))

    def __repr__(self):
        return "<Leaf: '{}'>".format(self.relative)
