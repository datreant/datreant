"""
Basic Container objects: the organizational units for :mod:`MDSynthesis`.

"""

import os
import yaml
import pdb
import shutil
import MDAnalysis

from Core import *

class Sim(ContainerCore):
    """The MDSynthesis Sim object is the base container for single simulations.

    The Sim object contains all the machinery required to handle trajectories
    and the data generated from them in an organized and object-oriented fashion.
    It is built to be used as a parent class, expanded upon for specific use
    cases as needed for specific simulation systems.

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

        ./MDSynthesis/Sim/name

    The object name can be specified as a keyword, or generated automatically.

    This directory contains a metadata file (Sim.yaml) with all the information
    needed by the object to find its trajectories and other generated data. To
    alter it, you need only open it in a text editor. This file is the same
    dictionary as in::

        s.metadata

    You can reload the metadata from the file with ``s.refresh()``. If you make
    changes to the metadata attribute interactively, you can write to the file
    using ``s.save()``.

    To regenerate an existing Sim object, give a directory that contains a Sim
    object metadata file (self.__class__.__name__ + ".yaml") instead of a topology::

        s = Sim('./MDSynthesis/Sim/name')

    The Sim object will be back as it was before.

    Data from Analysis objects are stored in the object directory. Having generated
    data from an Analysis called 'Foo', one would reload it with::

        s.load('Foo')

    and access it with::

        s.analysis['Foo']

    The data can be unloaded with::

        s.unload('Foo')

    This is beneficial if the data is rather large, freeing up memory. See the
    documentation for :class:`MDSynthesis.Operators.Analysis` for more details
    on how this scheme works.

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
            *selections*
                dictionary with selection strings to be applied to structure; 
                these will be used to make persistent AtomGroups; the dictionary
                may be recursive, i.e. items may be more dictionaries; the
                structure of this dictionary will be reflected in the Sim's
                selections attribute
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
        
        self.universe = None                     # universe 'dock'
        self.selections = dict()                 # AtomGroup selections
        self._cache = dict()                     # cache path storage
        self._uname = None                       # name of loaded universe

        if (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)

    def __repr__(self):
        if self._uname in self._cache:
            out = "{}: '{}' | universe (cached): '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)
        else:
            out = "{}: '{}' | universe: '{}'".format(self.__class__.__name__, self.metadata['name'], self._uname)

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
        self._init_database(database, locate=True)

        # record universe
        self.metadata['universes'] = dict()
        self.add('main', *args, **kwargs)

        # record selections
        self.metadata['selections'] = kwargs.pop('selections', dict())
            
        # finish up and save
        self.save()
        self._start_logger()

        # finally, attach universe to object
        if not detached:
            self.attach(uname)

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing Sim object.
        
        """
        attach = kwargs.pop('attach', None)
        load = kwargs.pop('load', None)

        basedir = os.path.abspath(args[0])
        self.metadata['basedir'] = basedir
        
        # get metadata (overwrites basedir metadata)
        self.refresh()
        
        # update location of object if changed
        self.metadata['basedir'] = basedir

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self.save()

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
    
        self.save()

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
        self.save()

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
            self._build_selections()
            self._logger.info("Atom selections generated.".format(self.metadata['name'], universe))
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
        loc = os.path.join(location, "{}.{}".format(self.metadata['uuid'], universe))

        i = 1
        location = "{}.{}".format(loc, i)
        while os.path.exists(location):
            i += 1
            location = "{}.{}".format(loc, i)

        self._cache[universe] = dict()
        self._cache[universe]['location'] = location
        self.util.makedirs(location)

        # build cached structure and trajectory filenames
        structure = self.metadata['universes'][universe]['structure']
        trajectory = self.metadata['universes'][universe]['trajectory']

        structure_c = os.path.join(location, os.path.basename(structure))
        trajectory_c = [ os.path.join(location, os.path.basename(x)) for x in trajectory ]

        # check before we accidentally overwrite valuable data
        if (structure_c == structure) or (trajectory_c == trajectory):
            self._logger.warning("Aborting cache; cache location same as storage!")
            return

        self._cache[universe]['structure'] = structure_c
        self._cache[universe]['trajectory'] = trajectory_c

        # copy to cache
        self._logger.info("Caching trajectory to {}\nThis may take some time...".format(location))
        shutil.copy2(structure, structure_c)
        for traj, traj_c in zip(trajectory, trajectory_c):
            shutil.copy2(traj, traj_c)

        self.universe = MDAnalysis.Universe(structure_c, *trajectory_c)
        self._logger.info("Universe '{}' now cached.".format(universe))
    
    def decache(self):
        """Remove cached trajectory from cache; re-attach universe to original.

        This operation deletes cached files, but ONLY cached files. It will not
        delete original trajectories.

        """
        if not (self._uname in self._cache):
            self._logger.info("Universe '{}' not cached.".format(universe))
            return

        self.attach(self._uname, force=True)

        structure = self.metadata['universes'][universe]['structure']
        trajectory = self.metadata['universes'][universe]['trajectory']
        
        structure_c = self._cache[universe]['structure']
        trajectory_c = self._cache[universe]['trajectory']
        location_c = self._cache[universe]['location']

        # final safety feature before delete
        if (structure_c == structure) or (trajectory_c == trajectory):
            self._logger.warning("Somehow cache is also storage location. Stopping before delete!")
        else:
            shutil.rmtree(location_c)

        del self._cache[universe]
        self._logger.info("Universe '{}' de-cached.".format(universe))

    def _build_selections(self):
        """Build selections attribute from selection strings stored in metadata.

        """
        def selection2atomGroup(selection):
            if isinstance(selection, dict):
                for key in selection:
                    selection[key] = selection2atomGroup(selection[key])
            elif isinstance(selection, basestring):
                agroup = self.universe.selectAtoms(selection)
                return agroup

        self.selections = selection2atomGroup(self.metadata['selections'])

    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default None
        """
        super(Sim, self)._build_metadata(**kwargs)

        # building core items
        uuid = self._generate_uuid()
        attributes = {'selections': kwargs.pop('selections', dict()),
                      }

        for key in attributes:
            if not key in self.metadata:
                self.metadata[key] = attributes[key]
        
class Group(ContainerCore):
    """Base class for a grouping of simulation objects.

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
        self.save()

    def remove(self, *args):
        """Remove a member from the Group.

        :Arguments:
            *args*
                uuid of member in self.members to be removed from Group
        """
        for uuid in args:
            self.metadata['members'].pop(uuid)
            self.members.pop(uuid)
        self.save()

    def _update_members(self):
        """Update member attributes.

        """
        for container in self.members:
            self.metadata['members'][container]['name'] = self.members[container].metadata['name']

        self.save()
            
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
        self.save()
        self._start_logger()
    
    def _regenerate(self, *args, **kwargs):
        """Re-generate existing object.
        
        """
        basedir = os.path.abspath(args[0])
        self.metadata['basedir'] = basedir
        
        # get metadata (overwrites basedir metadata)
        self.refresh()

        # update location of object if changed
        self.metadata['basedir'] = basedir

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self.save()

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

