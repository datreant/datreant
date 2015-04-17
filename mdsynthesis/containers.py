"""
Basic Container objects: the organizational units for :mod:`mdsynthesis`.

"""
import os
import sys
import shutil
from uuid import uuid4
import logging
from MDAnalysis import Universe

import core

class MultipleContainersError(Exception):
    pass

class NoContainersError(Exception):
    pass

class Container(object):
    """Core class for all Containers.

    """

    def __init__(self, container, location='.', coordinator=None,
            categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) generic Container
        object.

        :Required arguments:
            *container*
                if generating a new Container, the desired name to give it;
                if regenerating an existing Container, string giving the path
                to the directory containing the Container object's state file

        :Optional arguments when generating a new Container:
            *location*
                directory to place Container object; default is the current
                directory 
            *coordinator*
                directory of the Coordinator to associate with the Container;
                if the Coordinator does not exist, it is created; if ``None``,
                the Container will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give
                Containers distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing
                Container

        """
        if os.path.exists(container):
            self._regenerate('Container', container)
        else:
            self._generate('Container', container, location=location,
                    coordinator=coordinator, categories=categories, tags=tags)

    def _generate(self, containertype, container, location='.',
            coordinator=None, categories=None, tags=None):
        """Generate new generic Container object.
         
        """
        self._placeholders()

        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # generate state file
        #TODO: need try, except for case where Container already exists

        # TODO: need to raise exception where invalid characters used for directory
        os.makedirs(os.path.join(location, container))
        filename = core.filesystem.statefilename(containertype, str(uuid4()))

        statefile = os.path.join(location, container, filename)
        self._start_logger(containertype, container)

        self._containerfile = core.persistence.containerfile(statefile, self._logger,
                name=container, coordinator=coordinator, categories=categories,
                tags=tags)

    def _regenerate(self, containertype, container):
        """Re-generate existing Container object.
        
        """
        self._placeholders()

        # convenient to give only name of object (its directory name)
        if os.path.isdir(container):
            statefile = core.filesystem.glob_container(container)
            
            # if only one state file, load it; otherwise, complain loudly
            if len(statefile) == 1:
                self._containerfile = core.persistence.containerfile(statefile[0])
            elif len(statefile) == 0:
                raise NoContainersError('No Containers found in directory.')
            else:
                raise MultipleContainersError('Multiple Containers found in ' 
                        'directory. Give path to a specific state file.')

        # if a state file is given, try loading it
        elif os.path.exists(container):
            self._containerfile = core.persistence.containerfile(container)

        self._start_logger(containertype, self.name)
        self._containerfile._start_logger(self._logger)

        self._placeholders()

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        self._tags = None
        self._categories = None
        self._data = None

    def _start_logger(self, containertype, name, location=None, filehandler=False):
        """Start up the logger.

        :Arguments:
            *containertype*
                type of Container the logger is a part of; Sim or Group
            *name*
                name of Container
            *location*
                location of Container
            *filehandler*
                if True, write output to a logfile in the Container's main
                directory [``False``]

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(containertype, name))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)
    
            if filehandler:
                location = os.path.abspath(location)
                # file handler if desired; beware of problems with too many open files
                # when a large number of Containers are at play
                logfile = os.path.join(location, core.persistence.containerlog)
                fh = logging.FileHandler(logfile)
                ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
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
        if not os.path.exists(p):
            os.makedirs(p)

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

    @property
    def name(self):
        """The name of the Container.

        The name of a Container need not be unique with respect to other
        Containers, but is used as part of Container's displayed
        representation.
        
        """
        return os.path.basename(os.path.dirname(self._containerfile.filename))

    @name.setter
    def name(self, name):
        """The name of the Container.

        The name of a Container need not be unique with respect to other
        Containers, but is used as part of Container's displayed
        representation.
        
        """
        olddir = os.path.dirname(self._containerfile.filename)
        newdir = os.path.join(os.path.dirname(olddir), name)
        statefile = os.path.join(newdir,
                core.filesystem.statefilename(self.containertype, self.uuid))

        os.rename(olddir, newdir)
        self._regenerate(self.containertype, statefile)

    @property
    def containertype(self):
        """The type of the Container.
    
        """
        return os.path.basename(self._containerfile.filename).split('.')[0]

    @property
    def location(self):
        """The location of the Container.

        Setting the location to a new path physically moves the Container to
        the given location. This only works if the new location is an empty or
        nonexistent directory.
    
        """
        return os.path.dirname(self._containerfile.get_location())

    @location.setter
    def location(self, value):
        """Set location of Container. 
        
        Physically moves the Container to the given location.
        Only works if the new location is an empty or nonexistent
        directory.

        """
        self._makedirs(value)
        oldpath = self._containerfile.get_location()
        newpath = os.path.join(value, self.name)
        statefile = os.path.join(newpath,
                core.filesystem.statefilename(self.containertype, self.uuid))
        os.rename(oldpath, newpath)
        self._regenerate(self.containertype, statefile)
    
    @property
    def basedir(self):
        """Absolute path to the Container's base directory.

        This is a convenience property; the same result can be obtained by
        joining `:attr:location` and `:attr:name`.

        """
        return self._containerfile.get_location()

    @property
    def coordinator(self):
        """The location of the associated Coordinator.

        Change this to associate the Container with an existing
        or new Coordinator.
    
        """
        return self._containerfile.get_coordinator()

    #TODO: implement with Coordinator checking
    @coordinator.setter
    def coordinator(self, value):
        """Set location of Coordinator. 
        
        Setting this to ``None`` will dissociate the Container from any
        Coordinator. 
        
        """
        pass

    @property
    def tags(self):
        """The tags of the Container.

        Tags are user-added strings that can be used to and distinguish
        Containers from one another through Coordinator or Group queries.
        They can also be useful as flags for external code to determine
        how to handle the Container.
        
        """
        if not self._tags:
            self._tags = core.aggregators.Tags(self, self._containerfile, self._logger)
        return self._tags

    @property
    def categories(self):
        """The categories of the Container.
        
        Categories are user-added key-value pairs that can be used to and
        distinguish Containers from one another through Coordinator or Group
        queries. They can also be useful as flags for external code to
        determine how to handle the Container.

        """

        if not self._categories:
            self._categories = core.aggregators.Categories(self, self._containerfile, self._logger)
        return self._categories

    @property
    def data(self):
        """The data of the Container.

        Data are user-generated pandas objects (e.g. Series, DataFrames), numpy
        arrays, or any pickleable python object that are stored in the Container
        for easy recall later.  Each data instance is given its own directory
        in the Container's tree.
        
        """
        if not self._data:
            self._data = core.aggregators.Data(self, self._containerfile, self._logger)
        return self._data

    @property
    def uuid(self):
        """Get Container uuid.

        A Container's uuid is used by other Containers to identify it. The uuid
        is given in the Container's state file name for fast filesystem
        searching. For example, a Sim object with state file::

            'Sim.7dd9305a-d7d9-4a7b-b513-adf5f4205e09.h5'

        has uuid::

            '7dd9305a-d7d9-4a7b-b513-adf5f4205e09'

        Changing this string will alter the Container's uuid. This is not
        generally recommended.
    
        :Returns:
            *uuid*
                unique identifier string for this Container
        """
        return self._containerfile.filename.split('.')[1]

    def _new_uuid(self):
        """Generate new uuid for Container.

        *Warning*: A Container's uuid is used by Groups to identify whether or
        not it is a member. Any Groups a Container is a part of will cease
        to recognize it when changed.

        """
        new_id = str(uuid4())
        oldfile = self._containerfile.filename
        olddir = os.path.dirname(self._containerfile.filename)
        newfile = os.path.join(olddir, 
                core.filesystem.statefilename(self.containertype, uuid))
        os.rename(oldfile, newfile)
        self._regenerate(self.containertype, olddir)

class Sim(Container):
    """The Sim object is an interface to data for single simulations.
    
    """

    def __init__(self, sim, universe=None, uname='main', location='.',
                 coordinator=None, categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) Sim object.

        :Required arguments:
            *sim*
                if generating a new Sim, the desired name to give it;
                if regenerating an existing Sim, string giving the path
                to the directory containing the Sim object's state file

        :Optional arguments when generating a new Sim:
            *uname*
                desired name to associate with universe; this universe
                will be made the default (can always be changed later)
            *universe*
                arguments usually given to an MDAnalysis Universe
                that defines the topology and trajectory of the atoms
            *location*
                directory to place Sim object; default is the current directory
            *coordinator*
                directory of the Coordinator to associate with the Sim; if the
                Coordinator does not exist, it is created; if ``None``, the Sim
                will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give Sims
                distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing Sim

        """
        if os.path.exists(sim):
            self._regenerate('Sim', sim)
        else:
            self._generate('Sim', sim, universe=universe, uname=uname,
                    location=location, coordinator=coordinator,
                    categories=categories, tags=tags)

    def __repr__(self):
        if not self._uname:
            out = "<Sim: '{}'>".format(self.name)
        else:
            out = "<Sim: '{}' | active universe: '{}'>".format(self.name, self._uname)

        return out

    def __cmp__(self, other):
        if self.name < other.name:
            out = -1
        elif self.name == other.name:
            out = 0
        elif self.name > other.name:
            out = +1
        return out
    
    @property
    def universe(self):
        """The active universe of the Sim.

        Universes are interfaces to raw simulation data. The Sim can store
        multiple universe definitions corresponding to different versions
        of the same simulation output (e.g. post-processed trajectories derived
        from the same raw trajectory). The Sim has at most one universe
        definition that is "active" at one time, with stored selections for
        this universe directly available via ``Sim.selections``.

        To have more than one universe available as "active" at the same time,
        generate as many instances of the Sim object from the same statefile on
        disk as needed, and make a universe active for each one.
    
        """
        #TODO: include check for changes to universe definition, not just
        # definition absence
        if self._uname in self._containerfile.list_universes():
            return self._universe
        elif not self._universe:
            self.universes.activate()
            return self._universe
        else:
            self._universe = None
            self._logger.info('This universe is no longer defined. It has been detached.')

    @property
    def universes(self):
        """Manage the defined universes of the Sim.

        Universes are interfaces to raw simulation data. The Sim can store
        multiple universe definitions corresponding to different versions
        of the same simulation output (e.g. post-processed trajectories derived
        from the same raw trajectory). The Sim has at most one universe
        definition that is "active" at one time, with stored selections for
        this universe directly available via ``Sim.selections``.

        The Sim can also store a preference for a "default" universe, which is
        activated on a call to ``Sim.universe`` when no other universe is active.
        
        """
        if not self._universes:
            self._universes = core.aggregators.Universes(self, self._containerfile, self._logger)
        return self._universes

    @property
    def selections(self):
        """Stored atom selections for the active universe.

        Useful atom selections can be stored for the active universe and
        recalled later. Selections are stored separately for each defined
        universe, since the same selection may require a different selection
        string for different universes.

        """
        # attach default universe if not attached, and only give results if a
        # universe is present thereafter
        if self.universe:
            if not self._selections:
                self._selections = core.aggregators.Selections(self, self._containerfile, self._logger)
            return self._selections

    def _generate(self, containertype, container, universe=None, uname='main', location='.',
            coordinator=None, categories=None, tags=None):
        """Generate new Sim object.
         
        """
        super(Sim, self)._generate(containertype, container, location=location,
                coordinator=coordinator, categories=categories, tags=tags)

        # add universe
        if (uname and universe):
            self.universes.add(uname, *universe)
            self.universes.default(uname)

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        super(Sim, self)._placeholders()

        self._universes = None
        self._selections = None
        self._universe = None     # universe 'dock'
        self._uname = None        # attached universe name 

class Group(Container):
    """The Group object is a collection of Sims and Groups.

    """
    def __init__(self, group, members=None, location='.', coordinator=None, categories=None,
                 tags=None):
        """Generate a new or regenerate an existing (on disk) Group object.

        :Required Arguments:
            *group*
                if generating a new Group, the desired name to give it;
                if regenerating an existing Group, string giving the path
                to the directory containing the Group object's state file

        :Optional arguments when generating a new Group:
            *members*
                a list of Sims and/or Groups to immediately add as members
            *location*
                directory to place Group object; default is the current directory
            *coordinator*
                directory of the Coordinator to associate with this object; if the
                Coordinator does not exist, it is created; if ``None``, the Sim
                will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give
                Groups distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing Group

        """
        if os.path.exists(group):
            self._regenerate('Group', group)
        else:
            self._generate('Group', group, members=members, location=location,
                    coordinator=coordinator, categories=categories, tags=tags)

    def __repr__(self):
        members = self._containerfile.get_members_containertype()

        sims = members.count('Sim')
        groups = members.count('Group')

        out = "<Group: '{}'".format(self.name, len(members))
        if members:
            out = out +" | {} Members: ".format(len(members))
            if sims:
                out = out + "{} Sim".format(sims)
                if groups:
                    out = out + ", {} Group".format(groups)
            elif groups:
                out = out + "{} Group".format(groups)

        out = out + ">"

        return out

    @property
    def members(self):
        """The members of the Group.

        A Group is useful as an interface to collections of Containers, and
        they allow direct access to each member of that collection. Often
        a Group is used to store datasets derived from this collection as
        an aggregate.
        
        Queries can also be made on the Group's members to return a
        subselection of the members based on some search criteria. This can be
        useful to define new Groups from members of existing ones.
        
        """
        if not self._members:
            self._members = core.aggregators.Members(self, self._containerfile, self._logger)
        return self._members

    def _generate(self, containertype, container, members=None, location='.', coordinator=None,
                  categories=None, tags=None):
        """Generate new Group.
         
        """
        super(Group, self)._generate(containertype, container, location=location,
                coordinator=coordinator, categories=categories, tags=tags)

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
        self._cache = dict()    # member cache
