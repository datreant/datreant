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

        # quick fix for the case in which Database not found on generation
        if not 'database' in self.metadata:
            return None

    def __repr__(self):
        if self._uname in self._cache:
            out = "{}(Sim): '{}' | universe (cached): '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)
        else:
            out = "{}(Sim): '{}' | universe: '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)

        return out

    def info(self):
        """Output the current status of the Sim.

        """
        title = "{}: '{}'".format(self.__class__.__name__, self.metadata['name'])

        universes = "universes:\t"
        for universe in self.metadata['universes']:
            attached = cached = " "
            if universe in self._cache:
                cached = "(cached)"
            if self._uname == universe:
                attached = "*"

            universes = universes + "{} {} {}\n".format(universe, attached, cached)
            universes = universes + "\t\t"

        data = "data: "
        for datum in self.metadata['data']:
            loaded = " "
            if datum in self.data:
                loaded = "*"

            data = data + "\t\t{} {}\n".format(datum, loaded)

        out = "{}\n{}\n{}".format(title, universes, data)

        print out

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
    
    def add(self, name, *args, **kwargs):
        """Add a universe to the Sim.

        It is often useful to convert a raw MD trajectory into forms that are
        convenient for different types of analysis. Given this, the Sim object
        can associate with multiple trajectories, with the idea that all of
        these trajectories are drawn from the same raw simulation. The same Sim
        can thus be used to analyze every processed trajectory, keeping the
        data for a single simulation together.

        This method adds one universe at a time. 

        :Arguments:
            *name*
                identifier to use for the new universe
            *args*
                arguments normally given to ``MDAnalysis.Universe``

        :Keywords:
            **kwargs passed to ``MDAnalysis.Universe``

        """
        system = MDAnalysis.Universe(*args, **kwargs)
        self.metadata['universes'][name] = dict()

        self.metadata['universes'][name]['structure'] = os.path.abspath(system.filename)
        try:
            self.metadata['universes'][name]['trajectory'] = [ os.path.abspath(x) for x in system.trajectory.filenames ] 
        except AttributeError:
            self.metadata['universes'][name]['trajectory'] = [os.path.abspath(system.trajectory.filename)]
    
        self._save()

    def remove(self, *universe):
        """Remove a universe from Sim.

        This just removes the universe reference from the Sim's metadata, so
        it cannot be attached. This does not remove any other data, even if it
        was generated from the removed universe.

        :Arguments:
            *universe*
                identifier(s) of universe to remove
        """
        for n in universe:
            self.metadata['universes'].pop(n)
            if n == self._uname:
                self.detach()
        self._save()

    def attach(self, universe, **kwargs):
        """Attach universe.

        If another universe is already attached, it is detached first.
    
        :Arguments:
            *universe*
                universe to attach

        :Keywords:
            *force*
                if True, reattach universe even if already loaded; [``False``]
        """
        force = kwargs.pop('force', False)

        if (universe != self._uname) or (force == True):
            self._logger.info("Attaching '{}'...".format(universe))

            # attach cached universe if found to be cached; otherwise, attach original
            if universe in self._cache:
                self._logger.info("Found '{}' in cache; attaching cached version (decache to load original).".format(universe))
                structure = self._cache[universe]['structure']
                trajectory = self._cache[universe]['trajectory']
                self.universe = MDAnalysis.Universe(structure, *trajectory)
            else:
                structure = self.metadata['universes'][universe]['structure']
                trajectory = self.metadata['universes'][universe]['trajectory']
                self.universe = MDAnalysis.Universe(structure, *trajectory) 

            self._uname = universe
            self._logger.info("'{}' attached to universe '{}'.".format(self.metadata['name'], universe))
        else:
            self._logger.info("Skipping re-attach of {}...".format(i))

    def detach(self):
        """Detach universe.

        """
        self.universe = None
        self._logger.info("'{}' detached from universe'{}'".format(self.metadata['name'], self._uname))
        self._uname = None
    
    def cache(self, location):
        """Copy trajectory to a temporary location, and attach it as a universe.

        Useful if running analysis on a trajectory on a networked filesystem.
        This method will physically copy the structure and trajectory for the
        current universe to the specified *location*. Once transferred, the
        universe will be re-attached from the new location.

        :Arguments:
            *location*
                path to cache directory; will be made if it does not exist
        """
        if self._uname in self._cache:
            self._logger.warning("Aborting cache; universe already cached.")
            return
            
        # build and store location so we can delete it later
        location = os.path.abspath(location)
        loc = os.path.join(location, "{}.{}".format(self.metadata['uuid'], self._uname))

        i = 1
        location = "{}.{}".format(loc, i)
        while os.path.exists(location):
            i += 1
            location = "{}.{}".format(loc, i)

        self.util.makedirs(location)

        # build cached structure and trajectory filenames
        structure = self.metadata['universes'][self._uname]['structure']
        trajectory = self.metadata['universes'][self._uname]['trajectory']

        structure_c = os.path.join(location, os.path.basename(structure))
        trajectory_c = [ os.path.join(location, os.path.basename(x)) for x in trajectory ]

        # check before we accidentally overwrite valuable data
        if (structure_c == structure) or (trajectory_c == trajectory):
            self._logger.warning("Aborting cache; cache location same as storage!")
            return

        # copy to cache
        self._logger.info("Caching trajectory to {}\nThis may take some time...".format(location))
        shutil.copy2(structure, structure_c)
        for traj, traj_c in zip(trajectory, trajectory_c):
            shutil.copy2(traj, traj_c)

        self._cache[self._uname] = dict()
        self._cache[self._uname]['location'] = location
        self._cache[self._uname]['structure'] = structure_c
        self._cache[self._uname]['trajectory'] = trajectory_c

        self.universe = MDAnalysis.Universe(structure_c, *trajectory_c)
        self._logger.info("Universe '{}' now cached.".format(self._uname))
    
    def decache(self):
        """Remove cached trajectory from cache; re-attach universe to original.

        This operation deletes cached files, but ONLY cached files. It will not
        delete original trajectories.

        """
        if not (self._uname in self._cache):
            self._logger.info("Universe '{}' not cached.".format(universe))
            return

        structure = self.metadata['universes'][self._uname]['structure']
        trajectory = self.metadata['universes'][self._uname]['trajectory']
        
        structure_c = self._cache[self._uname]['structure']
        trajectory_c = self._cache[self._uname]['trajectory']
        location_c = self._cache[self._uname]['location']

        self._cache.pop(self._uname)
        self.attach(self._uname, force=True)

        # final safety feature before delete
        if (structure_c == structure) or (trajectory_c == trajectory):
            self._logger.warning("Somehow cache is also storage location. Stopping before delete!")
        else:
            shutil.rmtree(location_c)

        self._logger.info("Universe '{}' de-cached.".format(self._uname))
        
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

    def load_members(self, *args, **kwargs):
        """Load data instances for individual members of group.
    
        :Arguments:
            *args*
                datasets to load into each member

        :Keywords:
            *all*
                if True, unload all data instances [``False``]
        """
        for member in self.members.values():
            member.load(*args)

    def unload_members(self, *args):
        """Unload data instances from individual members of group.

        :Arguments:
            *args*
                datasets to unload from each member

        :Keywords:
            *all*
                if True, unload all data instances from each member [``False``]
        """
        for member in self.members.values():
            member.unload(*args)

    def add(self, *args): 
        """Add a member to the Group.
    
        :Arguments:
            *args*
                Sim-derived objects to add to Group
        """
        for container in args:
            self.metadata['members'].append({'name': system.metadata['name'],
                                             'type': system.metadata['type'],
                                             'basedir': system.metadata['basedir']
                                            })
            self.members.append(system)
        self._save()

    def remove(self, *args):
        """Remove a member from the Group.

        :Arguments:
            *args*
                uuid of member in self.members to be removed from Group
        """
        for uuid in args:
            self.metadata['members'].pop(uuid)
            self.members.pop(uuid)
        self._save()

    def _update_members(self):
        """Update member attributes.

        """
        for container in self.members:
            self.metadata['members'][container]['name'] = self.members[container].metadata['name']

        self._save()
            
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

    def _attach_members(self, **kwargs):
        """Attach member to Group object.
            
        Keyword arguments passed to Sim-derived object __init__().

        """
        self.members = dict()
        for key in self.metadata['members']:
            entry = self.metadata['members'][key]
            Simtype = self._simtype(entry['class'])
            self.members[key] = Simtype(entry['basedir'], **kwargs)
    
    def _simtype(self, typestring):
        """Return Sim or Sim-derived object based on type recorded in object
           metadata.

        :Arguments:
            *typestring*
                string pulled from Sim-derived object metadata signifying the
                object type (e.g. for a plain Sim object, this is ``Sim``)
        """
        objectout = None
        if typestring == 'Sim':
            objectout = Sim
            
        return objectout

    def _find_member(self, *uuid):
        """Find a member that has gone missing by consulting the Database.

        """
