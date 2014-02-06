"""
Basic Container objects: the organizational units for :mod:`MDSynthesis`.

"""

import os
import yaml
import pdb
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
                database does not exist, it is created
            *category*
                dictionary with user-defined keys and values; basically used to
                give Sims distinguishing characteristics
            *tag*
                list with user-defined values; like category, but useful for
                adding many distinguishing descriptors

        :Keywords always available:
            *detached*
                if True, Sim will load WITHOUT attaching trajectory; this is
                useful if only loadable analysis data are needed or
                trajectories are unavailable; default False

        """
        super(Sim, self).__init__()
        
        self.universe = dict()                  # universe 'modular dock'

        if (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)
    
    def _generate(self, *args, **kwargs):
        """Generate new Sim object.
         
        """
        system = MDAnalysis.Universe(*args, **kwargs)
        self._start_logger()

        # generate metadata items
        self._build_metadata(**kwargs)

        # find or generate database
        database = kwargs.pop('database', None)
        self._init_database(database, locate=True)

        # record universe
        self.metadata['universe']['main']['structure'] = os.path.abspath(system.filename)
        try:
            self.metadata['universe']['main']['trajectory'] = [ os.path.abspath(x) for x in system.trajectory.filenames ] 
        except AttributeError:
            self.metadata['universe']['main']['trajectory'] = [os.path.abspath(system.trajectory.filename)]
            
        # finish up and save
        self.save()
        self._start_logger()

        # finally, attach universe to object
        if not detached:
            self.attach('main')

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing Sim object.
        
        """
        detached = kwargs.pop('detached', False)
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

        # attach universe
        if not detached:
            self.attach('main')
    
    def attach(self, *args, **kwargs):
        """Attach universe.
    
        If 'all' is in argument list, every affiliated universe is loaded.

        :Keywords:
            *force*
                if True, reload data even if already loaded; default False
        """
        force = kwargs.pop('force', False)

        if 'all' in args:
            self._logger.info("Attaching all affiliated universes with '{}'...".format(self.metadata['name']))
            loadlist = self.metadata['universe']
        else:
            self._logger.info("Attaching selected universes to object '{}'...".format(self.metadata['name']))
            loadlist = args

        for i in loadlist:
            if (i not in self.universe) or (force == True):
                self._logger.info("Attaching {}...".format(i))
                structure = self.metadata['universe'][i]['structure']
                trajectory = [ x for x in self.metadata['universe'][i]['trajectory'] ]
                self.universe[i] = MDAnalysis.Universe(structure, *trajectory) 
            else:
                self._logger.info("Skipping re-attach of {}...".format(i))
        self._logger.info("Object '{}' attached to selected universes.".format(self.metadata['name']))

    def detach(self, *args, **kwargs):
        """Detach universe.

        If 'all' is in argument list, every loaded dataset is unloaded.

        :Arguments:
            *args*
                datasets to unload
        """
        if 'all' in args:
            self.universe.clear()
            self._logger.info("Object '{}' detached from all universes.".format(self.metadata['name']))
        else:
            self._logger.info("Detaching selected universes from object {}...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Detaching {}...".format(i))
                self.universe.pop(i, None)
            self._logger.info("Object '{}' detached from all selected universes.".format(self.metadata['name']))

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
            *tag*
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

    def load_members(self, *args):
        """Load data instances for individual members of group.
    
        If 'all' is in argument list, every available dataset is loaded for
        each member.

        :Arguments:
            *args*
                datasets to load into each member
        """
        for member in self.member:
            member.load(*args)

    def unload_members(self, *args):
        """Unload data instances from individual members of group.

        If 'all' is in argument list, every loaded dataset is unloaded from
        each member.

        :Arguments:
            *args*
                datasets to unload from each member
        """
        for member in self.member:
            member.unload(*args)

    def add(self, *args):
        """Add a member to the Group.

        :Arguments:
            *args*
                Sim-derived objects to add to Group
        """
        for system in args:
            self.metadata['member'].append({'name': system.metadata['name'],
                                             'type': system.metadata['type'],
                                             'basedir': system.metadata['basedir']
                                            })
            self.member.append(system)
        self.save()

    def remove(self, *args):
        """Remove a member from the Group.

        :Arguments:
            *args*
                index of member in self.members to be removed from Group
        """
        for index in args:
            self.metadata['member'].pop(index)
            self.member.pop(index)
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
        self.metadata['member'] = dict()
        for container in args:
            self.metadata['member'][container.metadata['uuid']] = {'name': container.metadata['name'],
                                                                    'class': container.metadata['class'],
                                                                    'basedir': container.metadata['basedir']
                                                                   }

        # attach members to object
        self.member = args

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

    def _attach_members(self, **kwargs):
        """Attach member to Group object.
            
        Keyword arguments passed to Sim-derived object __init__().

        """
        self.member = dict()
        for key in self.metadata['member']:
            entry = self.metadata['member'][key]
            Simtype = self._simtype(entry['class'])
            self.member[key] = Simtype(entry['basedir'], **kwargs)
    
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

