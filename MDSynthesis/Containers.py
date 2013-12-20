"""
Basic Container objects: the organizational units for :mod:`MDSynthesis`.

"""

import os
import yaml
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

        s.unload['Foo']

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

        :Keywords:
            *name*
                desired name for object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
            *projectdir*
                path to main project directory; defaults to current directory
            *pluck_segment*
                tuple with components of *trajpath* to leave out of final Sim
                object directory path, e.g. ('WORK/',)
            *naked*
                if True, Sim will load WITHOUT attaching trajectory or loading
                additional attributes; this is useful if only loadable analysis
                data are needed or trajectories are unavailable; default False
        """
        super(Sim, self).__init__()
        self.selections = dict()            # AtomGroups

        if (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)

    def _build_location(self, trajpath, *pluck_segment):
        """Build Sim object directory path from trajectory path.
    
        :Arguments:
            *trajpath*
                path to trajectory
            *pluck_segment*
                tuple with components of *trajpath* to leave out of final Sim
                object directory path, e.g. 'WORK/'
                
        """
        objectdir = '$PROJECT/MDSynthesis/{}'.format(self.__class__.__name__)

        # build path to container from trajpath
        p = os.path.abspath(trajpath)

        # pluck off trajectory filename from container path
        p = os.path.dirname(p)

        # pluck off projectdir part of path; replace with reference
        p = p.replace(self.metadata['projectdir'], objectdir)

        # subtract plucked segments from container path
        for seg in pluck_segment:
            seg = os.path.join(os.path.normpath(seg), '')
            p = p.replace(seg, '')

        # return final constructed path
        return p

    def _generate(self, *args, **kwargs):
        """Generate new Sim object.
         
        """
        system = MDAnalysis.Universe(*args, **kwargs)
        naked = kwargs.pop('naked', False)
        
        # set location of analysis structures
        projectdir = kwargs.pop('projectdir', None)
        if projectdir == None:
            # if no projectdir given, default to cwd
            self.metadata['projectdir'] = os.path.abspath('.')
        else:
            self.metadata['projectdir'] = os.path.abspath(projectdir)

        # process plucked segments
        pluck_segment = kwargs.pop('pluck_segment', ('',))
        if isinstance(pluck_segment, basestring):
            pluck_segment = [pluck_segment]
        else:
            pluck_segment = list(pluck_segment)
        self.metadata["basedir"] = self._build_location(system.trajectory.filename, *pluck_segment)
        
        self.metadata['metafile'] = '{}.yaml'.format(self.__class__.__name__)
        self.metadata['structure_file'] = self._abs2relpath(os.path.abspath(system.filename))

        # record trajectory file(s)
        try:
            self.metadata['trajectory_files'] = [ self._abs2relpath(os.path.abspath(x)) for x in system.trajectory.filenames ] 
        except AttributeError:
            self.metadata['trajectory_files'] = [self._abs2relpath(os.path.abspath(system.trajectory.filename))]

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self.save()

        # finally, attach universe to object
        if naked == False:
            self.universe = system
            self._build_attributes()

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing object.
        
        """
        naked = kwargs.pop('naked', False)
        basedir = os.path.abspath(args[0])
        metafile = os.path.join(basedir, '{}.yaml'.format(self.__class__.__name__))
        with open(metafile, 'r') as f:
            self.metadata = yaml.load(f)
        
        # update location of object if changed
        self._update_projectdir(basedir)
        self.metadata['basedir'] = self._abs2relpath(basedir)

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self.save()

        # attach universe
        if naked == False:
            self._attach_universe()
            self._build_attributes()
    
    def _attach_universe(self):
        """Attach universe, even if already attached.

        """
        structure = self._rel2abspath(self.metadata['structure_file'])
        trajectory = [ self._rel2abspath(x) for x in self.metadata['trajectory_files'] ]

        # attach universe
        self.universe = MDAnalysis.Universe(structure, *trajectory) 

    def _detach_universe(self):
        """Detach universe.

        """
        del self.universe

    def _build_metadata(self, **kwargs):
        """Build metadata. Runs on object generation. 
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
        """
        # fix name if object generated with no name or projectdir
        name = os.path.basename(os.path.dirname(self.metadata['trajectory_files'][0]))
        if name == '$PROJECT':
            name = self.__class__.__name__
        name = kwargs.pop('name', name)

        super(Sim, self)._build_metadata(name=name, **kwargs)

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
            *naked*
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
        for member in self.members:
            member.load(*args)

    def unload_members(self, *args):
        """Unload data instances from individual members of group.

        If 'all' is in argument list, every loaded dataset is unloaded from
        each member.

        :Arguments:
            *args*
                datasets to unload from each member
        """
        for member in self.members:
            member.unload(*args)

    def _generate(self, *args, **kwargs):
        """Generate new Group.
         
        """
        # get project directory from first object
        self.metadata['projectdir'] = args[0].metadata['projectdir']
        self.metadata['metafile'] = '{}.yaml'.format(self.__class__.__name__)

        # build list of Group members
        # will need to add unique hash references later
        self.metadata['members'] = list()
        for system in args:
            self.metadata['members'].append({'name': system.metadata['name'],
                                             'type': system.metadata['type'],
                                             'basedir': system.metadata['basedir']
                                            })

        self._build_metadata(**kwargs)
        self.metadata['basedir'] = '$PROJECT/MDSynthesis/{}/{}'.format(self.__class__.__name__, self.metadata['name'])

        # attach members to object
        self.members = args

        # finish up and save
        self._start_logger()
        self.save()
        self._build_attributes()

    def _regenerate(self, *args, **kwargs):
        """Re-generate existing object.
        
        """
        basedir = os.path.abspath(args[0])
        metafile = os.path.join(basedir, '{}.yaml'.format(self.__class__.__name__))
        with open(metafile, 'r') as f:
            self.metadata = yaml.load(f)
        
        # update location of object if changed
        self._update_projectdir(basedir)
        self.metadata['basedir'] = self._abs2relpath(basedir)

        # attach members to object
        self._attach_members(**kwargs)

        # finish up and save
        self._build_metadata(**kwargs)
        self._start_logger()
        self.save()
        self._build_attributes()

    def _attach_members(self, **kwargs):
        """Attach members to Group object.
            
        Keyword arguments passed to Sim-derived object __init__().

        """
        self.members = list()
        for entry in self.metadata['members']:
            Simtype = self._simtype(entry['type'])
            self.members.append(Simtype(self._rel2abspath(entry['basedir']),
                        **kwargs))
    
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
