"""
Basic Treant objects: the organizational units for :mod:`datreant`.

"""
import os
import sys
import shutil
from uuid import uuid4
import logging
import functools

import datreant
from .backends import statefiles
from . import limbs
from . import filesystem
from . import collections


class MultipleTreantsError(Exception):
    pass


class NoTreantsError(Exception):
    pass


def register(*treantclass):
    """Register a treant-derived class so datreant can handle it.

    In order for things like Bundle to know how to handle custom treants,
    they must be registered with the package. Give the class (not an instance)
    to this function to register it.

    :Arguments:
        *treantclass*
            treant-derived class to register; may enter more than one

    """
    for tc in treantclass:
        datreant._treants.update({tc._treanttype: tc})


@functools.total_ordering
class Treant(object):
    """Core class for all Treants.

    """
    # required components
    _treanttype = 'Treant'
    _backendclass = statefiles.TreantFile

    def __init__(self, treant, new=False, coordinator=None,
                 categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) Treant object.

        :Required arguments:
            *treant*
                base directory of a new or existing Treant; will regenerate
                a Treant if a state file is found, but will genereate a new
                one otherwise

                if multiple Treant state files are in the given directory,
                will raise :exception:`MultipleTreantsError`; specify
                the full path to the desired state file to regenerate the
                desired Treant in this case

                use the *new* keyword to force generation of a new Treant
                at the given path

        :Optional arguments:
            *new*
                generate a new Treant even if one already exists at the given
                location *treant*
            *categories*
                dictionary with user-defined keys and values; used to give
                Treants distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
        """
        if new:
            self._generate(treant, coordinator=coordinator,
                           categories=categories, tags=tags)
        else:
            try:
                self._regenerate(treant,
                                 coordinator=coordinator,
                                 categories=categories, tags=tags)
            except NoTreantsError:
                self._generate(treant,
                               coordinator=coordinator, categories=categories,
                               tags=tags)

        # any other init for this class of Treant
        self._init_hook()

    @classmethod
    def _attach_limb(cls, limb):
        """Attach a limb to the class, or to self if an instance.

        """
        # property definition
        def getter(self):
            if not hasattr(self, "_"+limb._name):
                setattr(self, "_"+limb._name, limb(self))
            return getattr(self, "_"+limb._name)

        # set the property
        setattr(cls, limb._name,
                property(getter, None, None, limb.__doc__))

    def __repr__(self):
        return "<Treant: '{}'>".format(self.name)

    def __getstate__(self):
        return self.filepath

    def __setstate__(self, state):
        self._regenerate(state)

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

    def __add__(a, b):
        """Addition of treants with collections or treants yields Bundle.

        """
        if (isinstance(a, (Treant, collections.CollectionBase)) and
           isinstance(b, (Treant, collections.CollectionBase))):
            return collections.Bundle(a, b)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def _generate(self, treant, coordinator=None, categories=None, tags=None):
        """Generate new Treant object.

        """

        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # build basedir; stop if we hit a permissions error
        try:
            self._makedirs(treant)
        except OSError as e:
            if e.errno == 13:
                raise OSError(13, "Permission denied; " +
                              "cannot create '{}'".format(treant))
            else:
                raise

        filename = filesystem.statefilename(self._treanttype, str(uuid4()))

        statefile = os.path.join(treant, filename)
        self._start_logger(self._treanttype, treant)

        # generate state file
        self._backend = datreant.backends.treantfile(
                statefile, self._logger, coordinator=coordinator,
                categories=categories, tags=tags)

    def _regenerate(self, treant, coordinator=None, categories=None,
                    tags=None):
        """Re-generate existing Treant object.

        """

        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # convenient to give only name of object (its directory name)
        if os.path.isdir(treant):
            statefile = filesystem.glob_treant(treant)

            # if only one state file, load it; otherwise, complain loudly
            if len(statefile) == 1:
                self._backend = datreant.backends.treantfile(
                        statefile[0], coordinator=coordinator,
                        categories=categories, tags=tags)
            elif len(statefile) == 0:
                raise NoTreantsError('No Treants found in directory.')
            else:
                raise MultipleTreantsError('Multiple Treants found in '
                                           'directory. Give path to a '
                                           'specific state file.')

        # if a state file is given, try loading it
        elif os.path.exists(treant):
            self._backend = datreant.backends.treantfile(
                    treant, coordinator=coordinator, categories=categories,
                    tags=tags)

        else:
            raise NoTreantsError('No Treants found in path.')

        self._start_logger(self._treanttype, self.name)
        self._backend._start_logger(self._logger)

    def _start_logger(self, treanttype, name, location=None,
                      filehandler=False):
        """Start up the logger.

        :Arguments:
            *treanttype*
                type of Treant the logger is a part of
            *name*
                name of Treant
            *location*
                location of Treant
            *filehandler*
                if True, write output to a logfile in the Treant's main
                directory [``False``]

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(treanttype, name))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)

            if filehandler:
                location = os.path.abspath(location)
                # file handler if desired; beware of problems with too many
                # open files when a large number of Treants are at play
                logfile = os.path.join(location, datreant.backends.treantlog)
                fh = logging.FileHandler(logfile)
                ff = logging.Formatter('%(asctime)s %(name)-12s '
                                       '%(levelname)-8s %(message)s')
                fh.setFormatter(ff)
                self._logger.addHandler(fh)

            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)

    def _makedirs(self, p):
        """Make directories and all parents necessary.

        :Arguments:
            *p*
                directory path to make
        """
        try:
            os.makedirs(p)
        except OSError as e:
            # let slide errors that include directories already existing, but
            # catch others
            if e.errno == 17:
                pass
            else:
                raise

    def _init_hook(self):
        """Perform any Treant-specific init tasks. For subclasses in other
        packages.

        """
        pass

    @property
    def name(self):
        """The name of the Treant.

        The name of a Treant need not be unique with respect to other
        Treants, but is used as part of Treant's displayed
        representation.

        """
        return os.path.basename(os.path.dirname(self._backend.filename))

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
        self._makedirs(value)
        oldpath = self._backend.get_location()
        newpath = os.path.join(value, self.name)
        statefile = os.path.join(newpath,
                                 filesystem.statefilename(
                                     self._treanttype, self.uuid))
        os.rename(oldpath, newpath)
        self._regenerate(statefile)

    @property
    def basedir(self):
        """Absolute path to the Treant's base directory.

        This is a convenience property; the same result can be obtained by
        joining :attr:`location` and :attr:`name`.

        """
        return self._backend.get_location()

    @property
    def filepath(self):
        """Absolute path to the Treant's state file.

        """
        return self._backend.filename

    def _new_uuid(self):
        """Generate new uuid for Treant.

        *Warning*: A Treant's uuid is used by Groups to identify whether or
        not it is a member. Any Groups a Treant is a part of will cease
        to recognize it when changed.

        """
        new_id = str(uuid4())
        oldfile = self._backend.filename
        olddir = os.path.dirname(self._backend.filename)
        newfile = os.path.join(olddir,
                               filesystem.statefilename(
                                   self._treanttype, uuid))
        os.rename(oldfile, newfile)
        self._regenerate(newfile)


class Group(Treant):
    """The Group object is a collection of Treants and Groups.

    """
    # required components
    _treanttype = 'Group'
    _backendclass = statefiles.GroupFile

    def __repr__(self):
        members = list(self._backend.get_members_treanttype())

        treants = members.count('Treant')
        groups = members.count('Group')

        out = "<Group: '{}'".format(self.name, len(members))
        if members:
            out = out + " | {} Members".format(len(members))
        out = out + ">"

        return out
