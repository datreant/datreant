"""
Basic Container objects: the organizational units for :mod:`MDSynthesis`.

"""

import os, sys
import shutil
import logging
from MDAnalysis import Universe
import Core

class _ContainerCore(object):
    """Core class for all Containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """
    def __init__(self):
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
                logfile = os.path.join(location, Core.Files.containerlog)
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

    @property
    def _uuid(self):
        """The uuid of the Container.

        This is automatically generated when the Container is first created,
        and is unchangeable. It is used by other MDSynthesis objects to
        uniquely identify the Container.
    
        """
        return self._containerfile.get_uuid()

    @property
    def name(self):
        """The name of the Container.

        The name of a Container is the only immutable element stored in the
        object's state file (aside from its unique id). It need not be unique
        with respect to other Containers, but is used as part of Container's
        displayed representation.
        
        """
        return self._containerfile.get_name()

    @property
    def _containertype(self):
        """The type of the Container; either Group or Sim.
    
        """
        return self._containerfile.get_containertype()

    @property
    def location(self):
        """The location of the Container.

        Setting the location to a new path physically moves the Container to
        the given location. This only works if the new location is an empty or
        nonexistent directory.
    
        """
        return self._containerfile.get_location()

    @location.setter
    def location(self, value):
        """Set location of Container. 
        
        Physically moves the Container to the given location.
        Only works if the new location is an empty or nonexistent
        directory.

        """
        self._makedirs(value)
        os.rename(self._containerfile.get_location(), value)
        self._regenerate(value)
    
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
            self._tags = Core.Aggregators.Tags(self, self._containerfile, self._logger)
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
            self._categories = Core.Aggregators.Categories(self, self._containerfile, self._logger)
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
            self._data = Core.Aggregators.Data(self, self._containerfile, self._logger)
        return self._data

class Sim(_ContainerCore):
    """The Sim object is an interface to data for single simulations.

    A Sim object contains all the machinery required to handle trajectories and
    the data generated from them in an organized and object-oriented fashion.

    To generate a Sim object from scratch, we need only give it a name. This
    will be used to distinguish the Sim from other Sims, though it need not be
    unique. We can also give it a topology and/or trajectory files as we would
    to an MDAnalysis Universe::
        
        s = Sim('fluffy', universe=[topology, trajectory])

    This will create a directory ``name`` that contains a single file (``Sim.h5``).
    That file is a persistent representation of the Sim on disk. We can access
    trajectory data by way of an MDAnalysis Universe::

        s.universe

    It can also store selections by giving the usual inputs to
    ``Universe.selectAtoms``::

        s.selections.add('backbone', ['name CA', 'name C', 'name O1', 'name O2'])

    And the AtomGroup can be conveniently obtained with::

        s.selections['backbone']

    The Sim can also store custom data structures. These can be pandas objects
    (e.g. Series, DataFrame, Panel), numpy arrays, or other python objects::

        a = np.random.randn(100, 100)
        s.data.add('randomdata', a)

    This can be recalled later with::

        s.data['randomdata']

    The real strength of the Sim object is how it stores its information. Generating
    an object from scratch stores the information needed to re-generate it in the
    filesystem. To generate another instance of the same Sim, simply give the directory
    where the state file lives::

        s2 = Sim('fluffy/')

    The Sim object will give access to the universe, stored selections, and stored data
    as before.

    """

    def __init__(self, sim, universe=None, uname='main', location='.',
                 coordinator=None, categories=None, tags=None, copy=None):
        """Generate or regenerate a Sim object.

        :Required arguments:
            *sim*
                if generating a new Sim, the desired name to give it;
                if regenerating an existing Sim, string giving the path
                to the directory containing the Sim object's state file

        :Optional arguments:
            *uname*
                desired name to associate with universe; this universe
                will be made the default (can always be changed later)
            *universe*
                arguments usually given to an MDAnalysis Universe
                that defines the topology and coordinates of the atoms
            *location*
                directory to place Sim object; default is current directory
            *coordinator*
                directory of the Coordinator to associate with this object; if the
                Coordinator does not exist, it is created [``None``] 
            *categories*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
            *copy* [TODO]
                if a Sim given, copy the contents into the new Sim;
                elements copied include tags, categories, universes,
                selections, and data 

        """
        super(Sim, self).__init__()

        self._universes = None
        self._selections = None
        self._universe = None     # universe 'dock'
        self._uname = None        # attached universe name 
        self._cache = dict()      # cache path storage

        if (os.path.isdir(sim)):
            # if directory string, load existing object
            self._regenerate(sim, universe=universe, uname=uname,
                    categories=categories, tags=tags, copy=copy)
        else:
            self._generate(sim, universe=universe, uname=uname,
                    location=location, coordinator=coordinator,
                    categories=categories, tags=tags, copy=copy)

    def __repr__(self):
        if not self._uname:
            out = "<Sim: '{}'>".format(self._containerfile.get_name())
        elif self._uname in self._cache:
            out = "<Sim: '{}' | active universe (cached): '{}'>".format(self._containerfile.get_name(), self._uname)
        else:
            out = "<Sim: '{}' | active universe: '{}'>".format(self._containerfile.get_name(), self._uname)

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
            self._universes = Core.Aggregators.Universes(self, self._containerfile, self._logger)
        return self._universes

    @property
    def selections(self):
        """Stored atom selections for the active universe.

        Useful atom selections can be stored for the active universe and
        recalled later. Selections are stored separately for each defined
        universe, since the same selection may require a different selection
        string for different universes.

        """
        # attach default universe if not attached
        self.universe
        if not self._selections:
            self._selections = Core.Aggregators.Selections(self, self._containerfile, self._logger)
        return self._selections

    def _generate(self, sim, universe=None, uname='main', location='.',
            coordinator=None, categories=None, tags=None, copy=None):
        """Generate new Sim object.
         
        """
        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # generate state file
        #TODO: need try, except for case where Sim already exists

        # name mangling to give a valid directory name
        # TODO: is this robust? What other characters are problematic?
        dirname = sim.replace('/', '_')
        os.makedirs(os.path.join(location, dirname))
        statefile = os.path.join(location, dirname, Core.Files.simfile)

        self._start_logger('Sim', sim)
        self._containerfile = Core.Files.SimFile(statefile, self._logger,
                name=sim, coordinator=coordinator, categories=categories,
                tags=tags)

        # add universe
        if (uname and universe):
            self.universes.add(uname, *universe)
            self.universes.default(uname)

    def _regenerate(self, sim, universe=None, uname='main', categories=None,
            tags=None, copy=None):
        """Re-generate existing Sim object.
        
        """
        # process keywords
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # load state file object
        statefile = os.path.join(sim, Core.Files.simfile)
        self._containerfile = Core.Files.SimFile(statefile,
                categories=categories, tags=tags)

        self._start_logger('Sim', self._containerfile.get_name())
        self._containerfile._start_logger(self._logger)

        # add universe
        if (uname and universe):
            self.universes.add(uname, *universe)

class Group(_ContainerCore):
    """The Group object is a collection of Sims and Groups.

    A Group object keeps track of any number of Sims and Groups added to it as
    members, and it can store datasets derived from these objects in the same
    way as Sims.

    To generate a Group object from scratch, we need only give it a name. We
    can also give any number of Sim and/or Groups as an argument::

        g = Group('gruffy', members=[s1, s2, s3, g4, g5])

    The Group can store custom data structures the same way as Sims. These can
    be pandas objects (e.g. Series, DataFrame, Panel), numpy arrays, or other
    python objects::

        a = np.random.randn(100, 100)
        g.data.add('randomdata', a)

    This can be recalled later with::

        g.data['randomdata']

    The real strength of the Group object is how it stores its information. Generating
    an object from scratch stores the information needed to re-generate it in the
    filesystem. To generate another instance of the same Group, simply give the directory
    where the state file lives::

        g2 = Group('gruffy/')

    The Group object will give access to its members and stored data as before.

    """
    def __init__(self, group, members=None, location='.', coordinator=None, categories=None,
                 tags=None, copy=None):
        """Generate or regenerate a Group object.

        :Required Arguments:
            *group*
                if generating a new Group, the desired name to give it;
                if regenerating an existing Group, string giving the path
                to the directory containing the Group object's state file

        :Optional arguments:
            *members*
                a list of Sims and/or Groups to immediately add as members
            *location*
                directory to place Group object; default is current directory
            *coordinator*
                directory of the Coordinator to associate with this object; if the
                Coordinator does not exist, it is created [``None``] 
            *categories*
                dictionary with user-defined keys and values; basically used to
                give Groups distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
            *copy* [TODO]
                if a Group given, copy the contents into the new Group;
                elements copied include tags, categories, members, and data

        """
        super(Group, self).__init__()
        self._members = None
        self._cache = dict()    # member cache

        if (os.path.isdir(group)):
            # if directory string, load existing object
            self._regenerate(group, members=members, categories=categories,
                    tags=tags, copy=copy)
        else:
            self._generate(group, members=members, location=location,
                    coordinator=coordinator, categories=categories, tags=tags,
                    copy=copy)

    def __repr__(self):
        members = self._containerfile.get_members_containertype()

        sims = members.count('Sim')
        groups = members.count('Group')

        out = "<Group: '{}' | {} Members: ".format(self._containerfile.get_name(), 
                                                len(members))
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
            self._members = Core.Aggregators.Members(self, self._containerfile, self._logger)
        return self._members

    def _generate(self, group, members=None, location='.', coordinator=None,
                  categories=None, tags=None, copy=None):
        """Generate new Group.
         
        """
        # process keywords
        if not members:
            members = list()
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # name mangling to give a valid directory name
        # TODO: is this robust? What other characters are problematic?
        dirname = group.replace('/', '_')
        os.makedirs(os.path.join(location, dirname))
        statefile = os.path.join(location, dirname, Core.Files.groupfile)

        self._start_logger('Group', group)
        self._containerfile = Core.Files.GroupFile(statefile, self._logger,
                name=group, coordinator=coordinator, categories=categories,
                tags=tags)

        # add members
        self.members.add(*members)
    
    def _regenerate(self, group, members=None, categories=None, tags=None,
            copy=None):
        """Re-generate existing object.
        
        """
        # process keywords
        if not members:
            members = list()
        if not categories:
            categories = dict()
        if not tags:
            tags = list()

        # load state file object
        statefile = os.path.join(group, Core.Files.groupfile)
        self._containerfile = Core.Files.GroupFile(statefile,
                categories=categories, tags=tags)

        self._start_logger('Group', self._containerfile.get_name())
        self._containerfile._start_logger(self._logger)

        # add members
        self.members.add(*members)
