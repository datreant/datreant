"""
Basic Treant objects: the organizational units for :mod:`datreant`.

"""
import os
import sys
import shutil
from uuid import uuid4
import logging
import functools

from datreant import aggregators
from datreant import filesystem
from datreant import persistence
from datreant import collections


class MultipleTreantsError(Exception):
    pass


class NoTreantsError(Exception):
    pass


@functools.total_ordering
class Treant(object):
    """Core class for all Treants.

    """
    _treanttype = 'Treant'

    def __init__(self, treant, location='.', coordinator=None,
                 categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) generic Treant
        object.

        :Required arguments:
            *treant*
                if generating a new Treant, the desired name to give it;
                if regenerating an existing Treant, string giving the path
                to the directory containing the Treant object's state file

        :Optional arguments when generating a new Treant:
            *location*
                directory to place Treant object; default is the current
                directory
            *coordinator*
                directory of the Coordinator to associate with the Treant;
                if the Coordinator does not exist, it is created; if ``None``,
                the Treant will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give
                Treants distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing
                Treant

        """
        if os.path.exists(os.path.join(location, treant)):
            self._regenerate(self._treanttype, treant,
                             coordinator=coordinator, categories=categories,
                             tags=tags)
        else:
            self._generate(self._treanttype, treant, location=location,
                           coordinator=coordinator, categories=categories,
                           tags=tags)

    def __repr__(self):
        return "<Treant: '{}'>".format(self.name)

    def __getstate__(self):
        return self.filepath

    def __setstate__(self, state):
        self._regenerate(self._treanttype, state)

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

    def __getitem__(self, handle):
        """Get dataset corresponding to given handle.

        If dataset doesn't exist, ``None`` is returned.

        :Arguments:
            *handle*
                name of data to retrieve

        :Returns:
            *data*
                stored data; ``None`` if nonexistent
        """
        return self.data.__getitem__(handle)

    def __setitem__(self, handle, data):
        """Set dataset corresponding to given handle.

        A data instance must be either a pandas Series, DataFrame, or Panel
        object. If dataset doesn't exist, it is added. If a dataset already
        exists for the given handle, it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data to store; must be a pandas Series, DataFrame, or Panel

        """
        self.data.__setitem__(handle, data)

    def __delitem__(self, handle):
        """Remove a dataset.

        Note: the directory containing the dataset file (``Data.h5``) will NOT
        be removed if it still contains file after the removal of the dataset
        file.

        :Arguments:
            *handle*
                name of dataset to delete

        """
        self.data.__delitem__(handle)

    def __add__(a, b):
        """Addition of treants with collections or treants yields Bundle.

        """
        if (isinstance(a, (Treant, collections._CollectionBase)) and
           isinstance(b, (Treant, collections._CollectionBase))):
            return collections.Bundle(a, b)
        else:
            raise TypeError("Operands must be Treant-derived or Bundles.")

    def _generate(self, treanttype, treant, location='.',
                  coordinator=None, categories=None, tags=None):
        """Generate new generic Treant object.

        """
        self._placeholders()

        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # generate state file
        # TODO: need try, except for case where Treant already exists

        # TODO: need to raise exception where invalid characters used for
        # directory
        self._makedirs(os.path.join(location, treant))
        filename = filesystem.statefilename(treanttype, str(uuid4()))

        statefile = os.path.join(location, treant, filename)
        self._start_logger(treanttype, treant)

        self._backend = persistence.treantfile(
                statefile, self._logger, name=treant,
                coordinator=coordinator, categories=categories, tags=tags)

    def _regenerate(self, treanttype, treant,
                    coordinator=None, categories=None, tags=None):
        """Re-generate existing Treant object.

        """
        self._placeholders()

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
                self._backend = persistence.treantfile(
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
            self._backend = persistence.treantfile(
                    treant, coordinator=coordinator, categories=categories,
                    tags=tags)

        self._start_logger(treanttype, self.name)
        self._backend._start_logger(self._logger)

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        self._tags = None
        self._categories = None
        self._data = None

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
                logfile = os.path.join(location, persistence.treantlog)
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
        except OSError:
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
                                     self.treanttype, self.uuid))

        os.rename(olddir, newdir)
        self._regenerate(self.treanttype, statefile)

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
                                     self.treanttype, self.uuid))
        os.rename(oldpath, newpath)
        self._regenerate(self.treanttype, statefile)

    @property
    def basedir(self):
        """Absolute path to the Treant's base directory.

        This is a convenience property; the same result can be obtained by
        joining `:attr:location` and `:attr:name`.

        """
        return self._backend.get_location()

    @property
    def filepath(self):
        """Absolute path to the Treant's state file.

        """
        return self._backend.filename

    @property
    def coordinators(self):
        """The locations of the associated Coordinators.

        Change this to associate the Treant with an existing
        or new Coordinator(s).

        """
        return self._backend.get_coordinator()

    # TODO: implement with Coordinator checking
    @coordinators.setter
    def coordinators(self, value):
        """Set locations of Coordinators.

        Setting this to ``None`` will dissociate the Treant from any
        Coordinators.

        """
        raise NotImplementedError("Coordinators are not yet implemented. This"
                                  " is a placeholder")

    @property
    def tags(self):
        """The tags of the Treant.

        Tags are user-added strings that can be used to and distinguish
        Treants from one another through Coordinator or Group queries.
        They can also be useful as flags for external code to determine
        how to handle the Treant.

        """
        if not self._tags:
            self._tags = aggregators.Tags(
                    self, self._backend, self._logger)
        return self._tags

    @property
    def categories(self):
        """The categories of the Treant.

        Categories are user-added key-value pairs that can be used to and
        distinguish Treants from one another through Coordinator or Group
        queries. They can also be useful as flags for external code to
        determine how to handle the Treant.

        """

        if not self._categories:
            self._categories = aggregators.Categories(
                    self, self._backend, self._logger)
        return self._categories

    @property
    def data(self):
        """The data of the Treant.

        Data are user-generated pandas objects (e.g. Series, DataFrames), numpy
        arrays, or any pickleable python object that are stored in the
        Treant for easy recall later.  Each data instance is given its own
        directory in the Treant's tree.

        """
        if not self._data:
            self._data = aggregators.Data(self, self._backend, self._logger)
        return self._data

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
        return self._backend.filename.split('.')[1]

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
                                   self.treanttype, uuid))
        os.rename(oldfile, newfile)
        self._regenerate(self.treanttype, olddir)


class Group(Treant):
    """The Group object is a collection of Treants and Groups.

    """
    _treanttype = 'Group'

    def __init__(self, group, members=None, location='.', coordinator=None,
                 categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) Group object.

        :Required Arguments:
            *group*
                if generating a new Group, the desired name to give it;
                if regenerating an existing Group, string giving the path
                to the directory containing the Group object's state file

        :Optional arguments when generating a new Group:
            *members*
                a list of Treants and/or Groups to immediately add as members
            *location*
                directory to place Group object; default is the current
                directory
            *coordinator*
                directory of the Coordinator to associate with this object; if
                the Coordinator does not exist, it is created; if ``None``, the
                Treant will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give
                Groups distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing
                Group

        """
        if os.path.exists(group):
            self._regenerate('Group', group)
        else:
            self._generate('Group', group, members=members, location=location,
                           coordinator=coordinator, categories=categories,
                           tags=tags)

    def __repr__(self):
        members = list(self._backend.get_members_treanttype())

        treants = members.count('Treant')
        groups = members.count('Group')

        # TODO: needs to work for any type of Treant subclasses
        out = "<Group: '{}'".format(self.name, len(members))
        if members:
            out = out + " | {} Members: ".format(len(members))
            if treants:
                if treants == 1:
                    out = out + "{} Treant".format(treants)
                elif treants > 1:
                    out = out + "{} Treants".format(treants)
                if groups:
                    out = out + ", "
            if groups:
                if groups == 1:
                    out = out + "{} Group".format(groups)
                elif groups > 1:
                    out = out + "{} Groups".format(groups)

        out = out + ">"

        return out

    @property
    def members(self):
        """The members of the Group.

        A Group is useful as an interface to collections of Treants, and
        they allow direct access to each member of that collection. Often
        a Group is used to store datasets derived from this collection as
        an aggregate.

        Queries can also be made on the Group's members to return a
        subselection of the members based on some search criteria. This can be
        useful to define new Groups from members of existing ones.

        """
        if not self._members:
            self._members = aggregators.Members(self, self._backend,
                                                self._logger)
        return self._members

    def _generate(self, treanttype, treant, members=None, location='.',
                  coordinator=None, categories=None, tags=None):
        """Generate new Group.

        """
        super(Group, self)._generate(treanttype, treant,
                                     location=location,
                                     coordinator=coordinator,
                                     categories=categories, tags=tags)

        # process keywords
        if not members:
            members = list()

        # add members
        self.members.add(*members)

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        super(Group, self)._placeholders()

        self._members = None
