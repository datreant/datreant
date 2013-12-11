"""
Useful container objects for :mod:`MDAnalysis` :class:`MDAnalysis.Universe`.
These are the organizational units for :mod:`MDSynthesis`.

"""

import os
import sys
import cPickle
import yaml
import logging
import MDAnalysis
import pdb

class Sim(object):
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
        """
        self.metadata = dict()              # information about object; defines base object
        self.selections = dict()            # AtomGroups
        self.analysis = dict()              # analysis data 'modular dock'

        if (os.path.isdir(args[0])):
        # if first arg is a directory string, load existing object
            self._regenerate(*args, **kwargs)
        else:
        # if a structure and trajectory(s) are given, begin building new object
            self._generate(*args, **kwargs)

    def load(self, *args):
        """Load data instances into object.

        If 'all' is in argument list, every available dataset is loaded.

        :Arguments:
            *args*
                datasets to load as attributes
        """
        if 'all' in args:
            self._logger.info("Loading all known data into object '{}'...".format(self.metadata['name']))
            for i in self.metadata['analysis_list']:
                self._logger.info("Loading {}...".format(i))
                with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self._logger.info("Object '{}' loaded with all known data.".format(self.metadata['name']))
        else:
            self._logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Loading {}...".format(i))
                with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self._logger.info("Object '{}' loaded with selected data.".format(self.metadata['name']))

    def unload(self, *args):
        """Unload data instances from object.

        If 'all' is in argument list, every loaded dataset is unloaded.

        :Arguments:
            *args*
                datasets to unload
        """
        if 'all' in args:
            self.analysis.clear()
            self._logger.info("Object '{}' unloaded of all data.".format(self.metadata['name']))
        else:
            self._logger.info("Unloading selected data from object {}...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Unloading {}...".format(i))
                self.analysis.pop(i, None)
            self._logger.info("Object '{}' unloaded of selected data.".format(self.metadata['name']))

    def save(self):
        """Save base object metadata.

        """
        basedir = self._rel2abspath(self.metadata["basedir"])
        self._makedirs(basedir)

        with open(os.path.join(basedir, self.metadata['metafile']), 'w') as f:
            yaml.dump(self.metadata, f)

    def refresh(self):
        """Reloads metadata from file.

        """
        basedir = self._rel2abspath(self.metadata['basedir'])
        metafile = os.path.join(basedir, self.metadata['metafile'])
        with open(metafile, 'r') as f:
            self.metadata = yaml.load(f)

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

    def _update_projectdir(self, basedir_abs):
        """Update projectdir based on given (and current) absolute path of basedir.

        """
        self.metadata['projectdir'] = basedir_abs.partition('/MDSynthesis')[0]

    def _abs2relpath(self, abspath):
        """Return path to file relative to project path.
        
        """
        return abspath.replace(self.metadata['projectdir'], '$PROJECT')

    def _rel2abspath(self, relpath):
        """Return realpath given a path relative to project directory. The
            opposite of _project_relpath.
        """
        return relpath.replace('$PROJECT', self.metadata['projectdir'])
        
    def _makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)

    def _generate(self, *args, **kwargs):
        """Generate new Sim object.
         
        """
        system = MDAnalysis.Universe(*args, **kwargs)
        
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

        # finally, attach universe to object
        self.universe = system

        # finish up and save
        self._build_metadata(**kwargs)
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
        structure = self._rel2abspath(self.metadata['structure_file'])
        trajectory = [ self._rel2abspath(x) for x in self.metadata['trajectory_files'] ]

        # attach universe
        self.universe = MDAnalysis.Universe(structure, *trajectory) 

        # finish up and save
        self._build_metadata(**kwargs)
        self.save()
        self._build_attributes()
    
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

        # building core items
        attributes = {'name': kwargs.pop('name', name),
                      'logfile': '{}.log'.format(self.__class__.__name__),
                      'analysis_list': [],
                      'type': self.__class__.__name__,
                      }

        for key in attributes.keys():
            if not key in self.metadata:
                self.metadata[key] = attributes[key]

    def _build_attributes(self):
        """Build attributes. Needed each time object is generated.

        """
        self._start_logger()

    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))
        self._logger.setLevel(logging.INFO)

        # file handler
        logfile = self._rel2abspath(os.path.join(self.metadata['basedir'], self.metadata['logfile']))
        fh = logging.FileHandler(logfile)
        ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        fh.setFormatter(ff)
        self._logger.addHandler(fh)

        # output handler
        ch = logging.StreamHandler(sys.stdout)
        cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        ch.setFormatter(cf)
        self._logger.addHandler(ch)
        
class SimSet(object):
    """Base class for a set of simulation objects.

    """
