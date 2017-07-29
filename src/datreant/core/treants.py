"""
Treants: the organizational units for :mod:`datreant`.

"""
import os
import functools

from .collections import Bundle
from .trees import Tree
from .names import treantdir_name
from .util import makedirs


@functools.total_ordering
class Treant(Tree):
    """The Treant: a discoverable Tree with metadata.

    `treant` should be the directory of a new or existing Treant.
    An existing Treant will be used if a `.datreant` directory is found inside.
    If no `.datreant` directory is found inside, a new Treant will be created.

    A Tree object may also be used in the same way as a directory string.

    Parameters
    ----------
    treant : str or Tree
        Base directory of a new or existing Treant; may also be a Tree object
    categories : dict
        dictionary with user-defined keys and values; used to give
        Treants distinguishing characteristics
    tags : list
        list with user-defined strings; like categories, but useful for
        adding many distinguishing descriptors
    """

    def _make_treantdir(self):
        abspath = super(Treant, self).__getattribute__('_path').absolute()
        treantdir = abspath / treantdir_name

        if not treantdir.exists():
            # build datreant dir; stop if we hit a permissions error
            try:
                makedirs(treantdir, exist_ok=True)
            except OSError as e:
                if e.errno == 13:
                    raise OSError(13, "Permission denied; " +
                                  "cannot create '{}'".format(treantdir))
                else:
                    raise

    def __init__(self, treant, categories=None, tags=None):
        # if given a Tree, get path out of it
        if isinstance(treant, Tree):
            treant = treant.abspath

        super(Treant, self).__init__(treant)

        # make treantdir
        self._make_treantdir()

        if tags:
            self.tags.add(tags)
        if categories:
            self.categories.add(categories)

    def attach(self, *limbname):
        """Attach limbs by name to this Treant.

        """
        for ln in limbname:
            try:
                limb = _TREELIMBS[ln]
            except KeyError:
                try:
                    limb = _LIMBS[ln]
                except KeyError:
                    raise KeyError("No such limb '{}'".format(ln))
            else:
                self._attach_limb(limb)

    def __repr__(self):
        return "<Treant: '{}'>".format(self.name)

    def __getstate__(self):
        return self.abspath

    def __setstate__(self, state):
        self.__init__(state)

    def __hash__(self):
        return hash(self.abspath)

    def __eq__(self, other):
        try:
            return self.abspath == other.abspath
        except AttributeError:
            return NotImplemented

    def __lt__(self, other):
        try:
            return self.name < other.name
        except AttributeError:
            return NotImplemented

    def __add__(self, other):
        """Addition of treants with collections or treants yields Bundle.

        """
        if isinstance(other, (Treant, Bundle)):
            limbs = self.limbs | other.limbs
            return Bundle(self, other, limbs=limbs)
        else:
            raise TypeError("Operands must be Treants or Bundles.")

    @property
    def name(self):
        """The name of the Treant.

        """
        return super(Treant, self).name

    @property
    def _treantdir(self):
        return os.path.join(self.abspath, treantdir_name)
