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
import datreant.persistence.pytables
import datreant.persistence.yaml
from datreant import aggregators
from datreant import filesystem
from datreant import collections


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
    _backends = {'pytables': ['.h5', datreant.persistence.pytables.TreantFile],
                 'yaml': ['.yml', datreant.persistence.yaml.TreantFile]}

    def __init__(self, treant, new=False, coordinator=None,
                 categories=None, tags=None, backend='pytables'):
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
            *backend*
                backend to use for a new state file; only 'pytables' available;
                not used for object regeneration

        """
        if new:
            self._generate(treant, coordinator=coordinator,
                           categories=categories, tags=tags, backend=backend)
        else:
            try:
                self._regenerate(treant,
                                 coordinator=coordinator,
                                 categories=categories, tags=tags)
            except NoTreantsError:
                self._generate(treant,
                               coordinator=coordinator, categories=categories,
                               tags=tags, backend=backend)

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

    def _generate(self, treant, coordinator=None, categories=None, tags=None,
                  backend='pytables'):
        """Generate new Treant object.

        """
        self._placeholders()

        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # generate state file

        # TODO: need to raise exception where invalid characters used for
        # directory
        self._makedirs(treant)
        filename = filesystem.statefilename(self._treanttype, str(uuid4()),
                                            self._backends[backend][0])

        statefile = os.path.join(treant, filename)
        self._start_logger(self._treanttype, treant)

        self._backend = datreant.persistence.treantfile(
                statefile, self._logger, coordinator=coordinator,
                categories=categories, tags=tags)

    def _regenerate(self, treant, coordinator=None, categories=None,
                    tags=None):
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
                self._backend = datreant.persistence.treantfile(
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
            self._backend = datreant.persistence.treantfile(
                    treant, coordinator=coordinator, categories=categories,
                    tags=tags)

        else:
            raise NoTreantsError('No Treants found in path.')

        self._start_logger(self._treanttype, self.name)
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
                logfile = os.path.join(location, datreant.persistence.treantlog)
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
                                     self._treanttype, self.uuid,
                                     self._backends[self.backend[0]][0]))

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
    def backend(self):
        """Get Treant backend type and class.

        A Treant stores its state in a persistent backend. The type of backend
        is chosen on creation of the Treant.

        :Returns:
            *backend*
                tuple giving the backend type and class
        """
        ext = os.path.splitext(self._backend.filename)[-1]
        for backend in self._backends:
            if ext == self._backends[backend][0]:
                out = (backend, self._backends[backend][1])

        return out

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
                                     self._treanttype, self.uuid,
                                     self._backends[self.backend[0]][0]))
        os.rename(oldpath, newpath)
        self._regenerate(statefile)

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
                                   self._treanttype, uuid,
                                   self._backends[self.backend[0]][0]))
        os.rename(oldfile, newfile)
        self._regenerate(newfile)


class Group(Treant):
    """The Group object is a collection of Treants and Groups.

    """
    # required components
    _treanttype = 'Group'
    _backends = {'pytables': ['.h5', datreant.persistence.pytables.GroupFile]}

    def __repr__(self):
        members = list(self._backend.get_members_treanttype())

        treants = members.count('Treant')
        groups = members.count('Group')

        out = "<Group: '{}'".format(self.name, len(members))
        if members:
            out = out + " | {} Members".format(len(members))
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

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        super(Group, self)._placeholders()

        self._members = None
