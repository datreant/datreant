"""
Treants: the organizational units for :mod:`datreant`.

"""
import os
import functools
import six
from uuid import uuid4
from pathlib2 import Path

from . import limbs
from . import filesystem
from .collections import Bundle
from .trees import Tree
from .util import makedirs

from .backends.statefiles import treantfile, TreantFile
from . import _TREANTS, _TREELIMBS, _LIMBS


class MultipleTreantsError(Exception):
    pass


class NoTreantsError(Exception):
    pass


class _Treantmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        treanttype = classdict['_treanttype']
        _TREANTS[treanttype] = cls


@functools.total_ordering
class Treant(six.with_metaclass(_Treantmeta, Tree)):
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

    Use the `new` keyword to force generation of a new Treant at the given
    path.

    Parameters
    ----------
    treant : str or Tree
        Base directory of a new or existing Treant; will regenerate
        a Treant if a state file is found, but will genereate a new
        one otherwise; may also be a Tree object
    new : bool
        Generate a new Treant even if one already exists at the given
        location
    categories : dict
        dictionary with user-defined keys and values; used to give
        Treants distinguishing characteristics
    tags : list
        list with user-defined values; like categories, but useful for
        adding many distinguishing descriptors
    """
    # required components
    _treanttype = 'Treant'
    _backendclass = TreantFile

    def __init__(self, treant, new=False, categories=None, tags=None):
        # if given a Tree, get path out of it
        if isinstance(treant, Tree):
            treant = treant.abspath

        if new:
            self._generate(treant, categories=categories, tags=tags)
        else:
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
        return "<{}: '{}'>".format(self._treanttype, self.name)

    def __getstate__(self):
        return self.filepath

    def __setstate__(self, state):
        self.__init__(state)

    def __hash__(self):
        return hash(self.uuid)

    def __eq__(self, other):
        try:
            return (self.name + self.uuid) == (other.name + other.uuid)
        except AttributeError:
            return NotImplemented

    def __lt__(self, other):
        try:
            return (self.name + self.uuid) < (other.name + other.uuid)
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

        filename = filesystem.statefilename(self._treanttype, str(uuid4()))

        statefile = os.path.join(treant, filename)

        # generate state file
        self._backend = treantfile(statefile)

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

        # convenient to give only name of object (its directory name)
        if os.path.isdir(treant):
            statefile = filesystem.glob_treant(treant)

            # if only one state file, load it; otherwise, complain loudly
            if len(statefile) == 1:
                self._backend = treantfile(statefile[0])
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

            elif len(statefile) == 0:
                raise NoTreantsError('No Treants found in directory.')
            else:
                raise MultipleTreantsError('Multiple Treants found in '
                                           'directory. Give path to a '
                                           'specific state file.')

        # if a state file is given, try loading it
        elif os.path.exists(treant):
            self._backend = treantfile(treant)
            # try to add categories, tags
            if categories or tags:
                # try to add categories, tags in one go
                try:
                    with self._write:
                        if categories:
                            self.categories.add(categories)
                        if tags:
                            self.tags.add(tags)
                except (OSError, IOError):
                    pass
        else:
            raise NoTreantsError('No Treants found in path.')

    @property
    def name(self):
        """The name of the Treant.

        The name of a Treant need not be unique with respect to other
        Treants, but is used as part of Treant's displayed
        representation.

        """
        return os.path.basename(self._backend.get_location())

    @name.setter
    def name(self, name):
        """The name of the Treant.

        The name of a Treant need not be unique with respect to other
        Treants, but is used as part of Treant's displayed
        representation.

        """
        olddir = os.path.dirname(self._backend.filename)
        newdir = os.path.join(os.path.dirname(olddir), name)
        statefile = os.path.join(newdir,
                                 filesystem.statefilename(
                                     self._treanttype, self.uuid))

        os.rename(olddir, newdir)
        self._regenerate(statefile)

    @property
    def uuid(self):
        """Get Treant uuid.

        A Treant's uuid is used by other Treants to identify it. The uuid
        is given in the Treant's state file name for fast filesystem
        searching. For example, a Treant with state file::

            'Treant.7dd9305a-d7d9-4a7b-b513-adf5f4205e09.h5'

        has uuid::

            '7dd9305a-d7d9-4a7b-b513-adf5f4205e09'

        Changing this string will alter the Treant's uuid. This is not
        generally recommended.

        :Returns:
            *uuid*
                unique identifier string for this Treant
        """
        return os.path.basename(self._backend.filename).split('.')[1]

    @property
    def treanttype(self):
        """The type of the Treant.

        """
        return os.path.basename(self._backend.filename).split('.')[0]

    @property
    def location(self):
        """The location of the Treant.

        Setting the location to a new path physically moves the Treant to
        the given location. This only works if the new location is an empty or
        nonexistent directory.

        """
        return os.path.dirname(self._backend.get_location())

    @location.setter
    def location(self, value):
        """Set location of Treant.

        Physically moves the Treant to the given location.
        Only works if the new location is an empty or nonexistent
        directory.

        """
        makedirs(value)
        oldpath = self._backend.get_location()
        newpath = os.path.join(value, self.name)
        statefile = os.path.join(newpath,
                                 filesystem.statefilename(
                                     self._treanttype, self.uuid))
        os.rename(oldpath, newpath)
        self._regenerate(statefile)

    @property
    def path(self):
        """Treant directory as a :class:`pathlib2.Path`.

        """
        return Path(self._backend.get_location())

    @property
    def filepath(self):
        """Absolute path to the Treant's state file.

        """
        return self._backend.filename

    @property
    def tree(self):
        """This Treant's directory as a Tree.

        """
        return Tree(self.abspath, limbs=self.limbs)

    @property
    def state(self):
        with self._read:
            state = self._state
        return state


class Group(Treant):
    """A Treant with a persistent Bundle of other Treants.

    `treant` should be a base directory of a new or existing Group. An
    existing Group will be regenerated if a state file is found.  If no state
    file is found, a new Group will be created.

    A Tree object may also be used in the same way as a directory string.

    If multiple Treant/Group state files are in the given directory, a
    :exc:`MultipleTreantsError` will be raised; specify the full path to the
    desired state file to regenerate the desired Group in this case. It is
    generally better to avoid having multiple state files in the same
    directory.

    Use the `new` keyword to force generation of a new Group at the given
    path.

    Parameters
    ----------
    treant : str or Tree
        Base directory of a new or existing Group; will regenerate
        a Group if a state file is found, but will genereate a new
        one otherwise; may also be a Tree object
    new : bool
        Generate a new Group even if one already exists at the given
        location
    categories : dict
        dictionary with user-defined keys and values; used to give
        Groups distinguishing characteristics
    tags : list
        list with user-defined values; like categories, but useful for
        adding many distinguishing descriptors
    """
    # required components
    _treanttype = 'Group'

    def __repr__(self):
        out = "<Group: '{}'".format(self.name)
        with self._read:
            try:
                n_mems = len(self._state['members'])
            except KeyError:
                n_mems = 0

        if n_mems:
            out = out + " | {} Members".format(n_mems)
        out = out + ">"

        return out
