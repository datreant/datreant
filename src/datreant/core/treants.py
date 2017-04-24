"""
Treants: the organizational units for :mod:`datreant`.

"""
import os
import functools
import six
from pathlib2 import Path

from .collections import Bundle
from .trees import Tree
from .names import treantdir_name, treantfile_name
from .util import makedirs

from .backends.statefiles import TreantFile
from . import _TREELIMBS, _LIMBS


class NoTreantsError(Exception):
    pass


@functools.total_ordering
class Treant(Tree):
    """The Treant: a Tree with a state file.

    `treant` should be a base directory of a new or existing Treant. An
    existing Treant will be regenerated if a state file is found.
    If no state file is found, a new Treant will be created.

    A Tree object may also be used in the same way as a directory string.

    If multiple Treant state files are in the given directory, a
    :exc:`MultipleTreantsError` will be raised; specify the full path to the
    desired state file to regenerate the desired Treant in this case. It is
    generally better to avoid having multiple state files in the same
    directory.

    Parameters
    ----------
    treant : str or Tree
        Base directory of a new or existing Treant; will regenerate
        a Treant if a state file is found, but will genereate a new
        one otherwise; may also be a Tree object
    categories : dict
        dictionary with user-defined keys and values; used to give
        Treants distinguishing characteristics
    tags : list
        list with user-defined values; like categories, but useful for
        adding many distinguishing descriptors
    """

    def __init__(self, treant, categories=None, tags=None):
        # if given a Tree, get path out of it
        if isinstance(treant, Tree):
            treant = treant.abspath

        try:
            self._regenerate(treant, categories=categories, tags=tags)
        except NoTreantsError:
            self._generate(treant, categories=categories, tags=tags)

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

    @property
    def _state(self):
        return self._backend._state

    @property
    def _read(self):
        return self._backend.read()

    @property
    def _write(self):
        return self._backend.write()

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

    def _generate(self, treant, categories=None, tags=None):
        """Generate new Treant object.

        """
        # build basedir; stop if we hit a permissions error
        try:
            makedirs(treant)
        except OSError as e:
            if e.errno == 13:
                raise OSError(13, "Permission denied; " +
                              "cannot create '{}'".format(treant))
            else:
                raise

        # build datreant dir
        treantdir = os.path.join(treant, treantdir_name)
        makedirs(treantdir)

        # generate state file
        statefile = os.path.join(treantdir, treantfile_name)
        self._backend = TreantFile(statefile)

        # add categories, tags in one go; doubles as file init so there's
        # something there
        if not tags:
            tags = []
        if not categories:
            categories = {}

        with self._write:
            self.categories.add(categories)
            self.tags.add(tags)

    def _regenerate(self, treant, categories=None, tags=None):
        """Re-generate existing Treant object.

        """

        # we give only name of object (its directory name)
        if os.path.isdir(treant):
            statefile = os.path.join(treant, treantdir_name, treantfile_name)

            # if only one state file, load it; otherwise, complain loudly
            if os.path.exists(statefile):
                self._backend = TreantFile(statefile)
                # try to add categories, tags in one go
                if categories or tags:
                    try:
                        with self._write:
                            if categories:
                                self.categories.add(categories)
                            if tags:
                                self.tags.add(tags)
                    except (OSError, IOError):
                        pass

            else:
                raise NoTreantsError('No Treants found in directory.')
        else:
            raise NoTreantsError('No Treants found in directory.')

    @property
    def name(self):
        """The name of the Treant.

        """
        return super(Treant, self).name

    @property
    def path(self):
        """Treant directory as a :class:`pathlib2.Path`.

        """
        return Path(self._backend.get_location()).absolute().parent

    @property
    def tree(self):
        """This Treant's directory as a Tree.

        """
        return Tree(self.abspath, limbs=self.limbs)

    @property
    def _treantdir(self):
        return os.path.join(self.abspath, treantdir_name)

    @property
    def _treantfile(self):
        return os.path.join(self.abspath, treantdir_name, treantfile_name)
