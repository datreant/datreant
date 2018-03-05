"""
Treants: the organizational units for :mod:`datreant`.

"""
import os
import functools
import six
from pathlib2 import Path

from .collections import Bundle
from .trees import Tree
from .names import TREANTDIR_NAME
from .util import makedirs
from .metadata import Tags, Categories
from .exceptions import NotATreantError


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

    def __init__(self, treant, categories=None, tags=None):
        # if given a Tree, get path out of it
        if isinstance(treant, Tree):
            treant = treant.abspath

        super(Treant, self).__init__(treant)

        # attach metadata objects
        self._tags = Tags(self)
        self._categories = Categories(self)

        # make treantdir
        self._make_treantdir()

        if tags:
            self.tags.add(tags)
        if categories:
            self.categories.add(categories)

    def _make_treantdir(self):
        abspath = super(Treant, self).__getattribute__('_path').absolute()
        treantdir = abspath / TREANTDIR_NAME

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
            return Bundle(self, other)
        else:
            raise TypeError("Operands must be Treants or Bundles.")

    @property
    def name(self):
        """The name of the Treant.

        """
        return super(Treant, self).name

    @property
    def _treantdir(self):
        return os.path.join(self.abspath, TREANTDIR_NAME)

    @property
    def tags(self):
        try:
            return self._tags
        except FileNotFoundError:
            raise NotATreantError("This Treant no longer has a '.datreant'"
                                  " directory.")

    @tags.setter
    def tags(self, value):
        try:
            if isinstance(value, (Tags, list, set)):
                val = list(value)
                self.tags.clear()
                self.tags.add(val)
            else:
                raise TypeError("Can only set with tags, a list, or set")
        except FileNotFoundError:
            raise NotATreantError("This Treant no longer has a '.datreant'"
                                  " directory.")

    @property
    def categories(self):
        try:
            return self._categories
        except FileNotFoundError:
            raise NotATreantError("This Treant no longer has a '.datreant'"
                                  " directory.")

    @categories.setter
    def categories(self, value):
        try:
            if isinstance(value, (Categories, dict)):
                val = dict(value)
                self.categories.clear()
                self.categories.add(val)
            else:
                raise TypeError("Can only set with categories or dict")
        except FileNotFoundError:
            raise NotATreantError("This Treant no longer has a '.datreant'"
                                  " directory.")
