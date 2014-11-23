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
    def _start_logger(self, containertype, name, location):
        """Start up the logger.

        :Arguments:
            *containertype*
                type of Container the logger is a part of; Sim or Group
            *name*
                name of Container
            *location*
                location of Container

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(containertype, name))

        location = os.path.abspath(location)
        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)
    
            # file handler
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
    
        """
        return self._containerfile.get_uuid()

    @property
    def name(self):
        """The name of the Container.
        
        """
        return self._containerfile.get_name()

    @name.setter
    def name(self, value):
        """Set name of Container.

        """
        self._containerfile.update_name(value)

    @property
    def _containertype(self):
        """The type of the Container; either Group or Sim.
    
        """
        return self._containerfile.get_containertype()

    @property
    def location(self):
        """The location of the Container.
    
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

#TODO: include in documentation fgetter details
class Sim(_ContainerCore):
    """The MDSynthesis Sim object is the base container for single simulations.

    The Sim object contains all the machinery required to handle trajectories
    and the data generated from them in an organized and object-oriented fashion.

    To generate a Sim object from scratch, provide a topology and a trajectory
    in the same way you would for a Universe (:class:`MDAnalysis.Universe`). 

    For example, as with a Universe::

       s = Sim(topology, trajectory)          # read system from file(s)
       s = Sim(pdbfile)                       # read atoms and coordinates from PDB or GRO
       s = Sim(topology, [traj1, traj2, ...]) # read from a list of trajectories
       s = Sim(topology, traj1, traj2, ...)   # read from multiple trajectories

    The real strength of the Sim object is how it stores its information. Generating
    an object from scratch stores the information needed to re-generate it in the
    filesystem. By default, this is the current working directory::

        ./Sim

    This directory contains a state file with all the information needed by the
    object to find its trajectories and other generated data.

    To regenerate an existing Sim object, give a directory that contains a Sim
    object state file instead of a topology::

        s = Sim('path/to/sim/directory')

    The Sim object will be back as it was before.

    """

    def __init__(self, *args, **kwargs):
        """Generate or regenerate a Sim object.

        :Arguments:
              *args*
                either a string giving the path to a directory with a Sim object
                state file, or the arguments normally given to an MDAnalysis
                Universe

        :Keywords available on object generation:
            *location*
                directory to place Sim object; default is current directory
            *name*
                desired name for Sim; can be used to distinguish Sims, but need
                not be unique; default is ``Sim``
            *coordinator*
                directory of the Coordinator to associate with this object; if the
                Coordinator does not exist, it is created [``None``] 
            *universe*
                desired name to associate with first universe [``main``]
            *categories*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
            *detached*
                if True, Sim will load WITHOUT attaching Universe; default
                False 
            *empty*
                if True, initialize Sim without Universe definition;
                no arguments required; default False

        :Keywords available on object re-generation:
            *attach*
                name of universe to attach

        """

        self._universe = None     # universe 'dock'
        self._uname = None        # attached universe name 
        self._cache = dict()      # cache path storage

        if kwargs.pop('empty', False):
        # if no universe desired, skip checks for arguments
            self._generate(*args, empty=True, **kwargs)
        elif (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)

    def __repr__(self):
        if self._uname in self._cache:
            out = "Sim: '{}' | universe (cached): '{}'".format(self._containerfile.get_name(), self._uname)
        else:
            out = "Sim: '{}' | universe: '{}'".format(self._containerfile.get_name(), self._uname)

        return out

    @property
    def universe(self):
        """The active Universe of the Sim.
    
        """
        if self._uname in self._containerfile.list_universes():
            return self._universe
        elif not self._universe:
            self._logger.info('No Universe attached.')
        else:
            self.detach()
            self._logger.info('This Universe is no longer defined. It has been detached')

    #TODO: add explicit args, kwargs
    def _generate(self, *args, **kwargs):
        """Generate new Sim object.
         
        """
        # process keywords
        location = kwargs.pop('location', '.')
        name = kwargs.pop('name', 'Sim')
        coordinator = kwargs.pop('coordinator', None)
        categories = kwargs.pop('categories', dict())
        tags = kwargs.pop('tags', list())
        universe = kwargs.pop('universe', 'main')
        detached = kwargs.pop('detached', False)
        empty = kwargs.pop('empty', False)

        # generate state file
        #TODO: need try, except for case where Sim already exists

        # name mangling to give a valid directory name
        # TODO: is this robust? What other characters are problematic?
        dirname = name.replace('/', '_')
        os.makedirs(os.path.join(location, dirname))
        statefile = os.path.join(location, dirname, Core.Files.simfile)

        self._start_logger('Sim', name, os.path.join(location, dirname))
        self._containerfile = Core.Files.SimFile(statefile, self._logger,
                                                 name=name,
                                                 coordinator=coordinator,
                                                 categories=categories,
                                                 tags=tags)

        # attach aggregators
        self._init_aggregators()

        # add universe
        if not empty:
            self.universes.add(universe, args[0], *args[1:])
            if not detached:
                self.attach(universe)

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing Sim object.
        
        """
        attach = kwargs.pop('attach', None)

        #TODO: load logger first!

        # load state file object
        statefile = os.path.join(args[0], Core.Files.simfile)
        self._containerfile = Core.Files.SimFile(statefile)
        self._start_logger('Sim', self._containerfile.get_name(), args[0])
        self._containerfile._start_logger(self._logger)

        # attach aggregators
        self._init_aggregators()
    
        if attach:
            self.attach(attach)

    def _init_aggregators(self):
        """Initialize and attach aggregators.

        """
        #TODO: INFO needs work
        #self.info = Core.Aggregators.SimInfo(self, self._containerfile, self._logger)
        self.tags = Core.Aggregators.Tags(self, self._containerfile, self._logger)
        self.categories = Core.Aggregators.Categories(self, self._containerfile, self._logger)
        self.universes = Core.Aggregators.Universes(self, self._containerfile, self._logger)
        self.selections = Core.Aggregators.Selections(self, self._containerfile, self._logger)
        self.data = Core.Aggregators.Data(self, self._containerfile, self._logger)

    def attach(self, handle):
        """Attach the given universe.
        
        Only one universe definition can be attached to a Sim at one time. The
        attached universe can be accessed from ``Sim.universe``.
    
        :Arguments:
            *handle*
                given name for selecting the universe
        """
        udef = self._containerfile.get_universe(handle)
        self._universe = Universe(udef[0], *udef[1])
        self._uname = handle
    
    def detach(self):
        """Detach from universe.

        Deletes the currently attached Universe object from ``Sim.universe``.

        """
        self._universe = None
        self._uname = None

class Group(_ContainerCore):
    """

    """
    def __init__(self, *args, **kwargs):
        """Generate or regenerate a Group object.

        :Arguments:
              *args*
                either a string giving the path to a directory with a Group object
                metadata file, or any number of re-generated Sim-derived
                objects 

        :Keywords:
            *name*
                desired name for object, used for logging and referring to
                object in some analyses; default is class name
            *database*
                directory of the database to associate with this object; if the
                database does not exist, it is created; if none is specified, a
                database is created in the current directory
            *categories*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
            *detached*
                if True, members will load WITHOUT attaching trajectories or
                loading additional attributes; this is useful if only loadable
                analysis data are needed or trajectories are unavailable;
                default False
                
        """

        if isinstance(args[0], basestring):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a number of Sim-derived objects are given, build a new group 
            self._generate(*args, **kwargs)

    def _generate(self, *args, **kwargs):
        """Generate new Group.
         
        """
        # generate metadata items
        self._build_metadata(**kwargs)
        self._start_logger()

        # find or generate database
        database = kwargs.pop('database', None)
        self._init_database(database, locate=True)

        # build list of Group members
        self.metadata['members'] = dict()
        for container in args:
            self.metadata['members'][container.metadata['uuid']] = {'name': container.metadata['name'],
                                                                    'class': container.metadata['class'],
                                                                    'basedir': container.metadata['basedir']
                                                                   }

        # attach members to object
        self.members = args

        # finish up and save
        self._save()
        self._start_logger()
    
    def _regenerate(self, *args, **kwargs):
        """Re-generate existing object.
        
        """
        basedir = os.path.abspath(args[0])
        self.metadata['basedir'] = basedir
        
        # get metadata (overwrites basedir metadata)
        self._refresh()

        # update location of object if changed
        self.metadata['basedir'] = basedir

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self._save()

        # attach members to object
        self._attach_members(**kwargs)
        self._update_members()

