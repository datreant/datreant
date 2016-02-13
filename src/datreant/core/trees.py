import os
from six import string_types
import scandir
from pathlib import Path

from asciitree import LeftAligned

from .util import makedirs

class TreeMixin(object):
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
        fullpath = os.path.join(self.path, path)

        if os.path.isdir(fullpath) or fullpath.endswith('/'):
            return Tree(fullpath)
        else:
            return Leaf(fullpath)

    @property
    def leaves(self):
        """A list of the file names in the directory.

        """
        for root, dirs, files in scandir.walk(self.path):
            out = files
            break

        out.sort()
        return out

    @property
    def trees(self):
        """A list of the directories in the directory.

        """
        for root, dirs, files in scandir.walk(self.path):
            out = dirs 
            break

        out.sort()
        return out

    @property
    def draw(self):
        """Print an asciified visual of the tree.

        """
        tree = {}
        rootdir = self.path.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in os.walk(rootdir):
            folders = ["{}/".format(x) for x in path[start:].split(os.sep)]
            subdir = dict.fromkeys(files, {})
            parent = reduce(dict.get, folders[:-1], tree)
            parent[folders[-1]] = subdir
        
        tr = LeftAligned()
        print tr(tree)


class BrushMixin(object):
    def __str__(self):
        return str(self.path)

    @property
    def exists(self):
        return self.path.exists

    @property
    def path(self):
        return self._path

    @property
    def abspath(self):
        return str(self.path.absolute())

    @property
    def relpath(self):
        return str(self.path.relative_to(os.path.abspath('.')))


class Tree(TreeMixin, BrushMixin):
    """A directory.

    """

    def __init__(self, dirpath):
        makedirs(dirpath)
        self._path = Path(os.path.abspath(dirpath))

    def __repr__(self):
        return "<Tree: '{}'>".format(self.relative)


class Leaf(BrushMixin):
    """A file in the filesystem.

    """

    def __init__(self, filepath):
        makedirs(os.path.dirname(filepath))
        self._path = Path(os.path.abspath(filepath))

    def __repr__(self):
        return "<Leaf: '{}'>".format(self.relative)
