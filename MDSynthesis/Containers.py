"""
Basic Container objects: the organizational units for :mod:`MDSynthesis`.

"""

import os, sys
import shutil
import logging
import MDAnalysis
import Core

class _ContainerCore(Core.Workers.ObjectCore):
    """Core class for all Containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """

    def __init__(self):
        """
        
        """
        super(ContainerCore, self).__init__()

        #TODO: update location state

    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)
    
            # file handler
            if 'basedir' in self.metadata:
                logfile = os.path.join(self.metadata['basedir'], self._containerlog)
                fh = logging.FileHandler(logfile)
                ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
                fh.setFormatter(ff)
                self._logger.addHandler(fh)
    
            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)
    
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
                metadata file, or the arguments normally given to an MDAnalysis
                Universe

        :Keywords used on object generation:
            *name*
                desired name for object, used for logging and referring to
                object in some analyses; default will be the object's randomly
                selected UUID
            *database*
                directory of the database to associate with this object; if the
                database does not exist, it is created; if none given, the
                directory tree will be searched upward until one or none is found
            *universe*
                desired name to associate with first universe [``main``]
            *category*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tags*
                list with user-defined values; like category, but useful for
                adding many distinguishing descriptors
            *detached*
                if True, Sim will load WITHOUT attaching trajectory; this is
                useful if only loadable analysis data are needed or
                trajectories are unavailable; default False

        :Keywords used on object re-generation:
            *attach*
                name of universe to attach
            *load*
                list of data elements to load

        :Keywords always available:

        """
        super(Sim, self).__init__()

        self._containerfile = Core.Files.SimFile(Core.simfile)
        
        self.universe = None      # universe 'dock'
        self._uname = None        # attached universe name 
        self.selections = Core.Aggregators.Selections(self, self._containerfile, self._logger)
        self._cache = dict()      # cache path storage

        if (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)

    def __repr__(self):
        if self._uname in self._cache:
            out = "{}(Sim): '{}' | universe (cached): '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)
        else:
            out = "{}(Sim): '{}' | universe: '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)

        return out

    def _generate(self, *args, **kwargs):
        """Generate new Sim object.
         
        """
        detached = kwargs.pop('detached', False)
        uname = kwargs.pop('universe', 'main')
        system = MDAnalysis.Universe(*args, **kwargs)

        # generate metadata items
        self._build_metadata(**kwargs)
        self._start_logger()

        # find or generate database
        database = kwargs.pop('database', None)
        success = self._init_database(database, locate=True)

        if success:
            # record universe
            self.metadata['universes'] = dict()
            self.add('main', *args, **kwargs)
                
            # finish up and save
            self._save()
            self._start_logger()

            # finally, attach universe to object
            if not detached:
                self.attach(uname)
        else:
            self._logger.info("Terminating generation. Specify a database location.")

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing Sim object.
        
        """
        attach = kwargs.pop('attach', None)
        load = kwargs.pop('load', None)

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

        # attach universe, if desired
        if attach:
            self.attach(attach)
    
class Group(_ContainerCore):
    """A grouping of Sim objects.

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
            *category*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tags*
                list with user-defined values; like category, but useful for
                adding many distinguishing descriptors
            *detached*
                if True, members will load WITHOUT attaching trajectories or
                loading additional attributes; this is useful if only loadable
                analysis data are needed or trajectories are unavailable;
                default False
                
        """
        super(Group, self).__init__()

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

