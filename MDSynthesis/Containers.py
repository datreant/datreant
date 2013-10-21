"""
Useful container objects for :mod:`MDAnalysis` :class:`MDAnalysis.core.AtomGroup.Universe`.
These are the organizational units for :mod:`MDSynthesis`.

"""

import os
import sys
import cPickle
import yaml
import logging
import MDAnalysis

class Sim(object):
    """Base class for simulation objects.

    """

    def __init__(self, system, **kwargs):
        """Initialize a Sim object.

        :Arguments:
            *system*
                universe (:class:`MDAnalysis.Universe` object) that contains a
                trajectory; can also be a directory that contains a Sim object
                metadata file

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
            *location*
                location for analysis data storage; default is analysis
                directory corresponding to trajectory directory
                if *system* is a base object directory, location is set
                to that directory
            *projectdir*
                path to main project directory; required if no *location* given

        """
        self.metadata = dict()              # information about object; defines base object
        self.selections = dict()            # AtomGroups
        self.analysis = dict()              # analysis data 'modular dock'

        # if system is a directory string, load existing base object
        if isinstance(system, basestring):
            self.metadata["basedir"] = os.path.abspath(system)
            self.metadata["metafile"] = os.path.join(self.metadata["basedir"], '{}.yaml'.format(self.__name__))
            self._load_base()
        # if system is a universe, begin building new base object
        elif isinstance(system, MDAnalysis.core.AtomGroup.Universe):
            # set location of analysis structures
            location = kwargs.pop('location', None)
            if location == None:
                try:
                    projectdir = kwargs.pop('projectdir')
                except KeyError:
                    print "Cannot construct {} object without projectdir. See documentation for details.".format(self.__name__)
                analysisdir = os.path.join(projectdir, 'MDSynthesis/{}'.format(self.__name__))
                self.metadata["basedir"] = self._location(system.trajectory.filename, analysisdir)
            else:
                self.metadata["basedir"] = os.path.abspath(location)
            self.metadata["metafile"] = os.path.join(self.metadata["basedir"], '{}.yaml'.format(self.__name__))
            self.metadata['structure_file'] = os.path.abspath(system.filename) 
            self.metadata['trajectory_file'] = os.path.abspath(system.trajectory.filename)
            self.universe = system

        # finish up and save
        self._build_metadata(**kwargs)
        self.save()

    def load(self, *args):
        """Load data instances into object.

        If 'all' is in argument list, every available dataset is loaded.

        :Arguments:
            *args*
                datasets to load as attributes
        """
        if 'all' in args:
            self.logger.info("Loading all known data into object '{}'...".format(self.metadata['name']))
            for i in self.metadata['analysis_list']:
                self.logger.info("Loading {}...".format(i))
                with open(os.path.join(self.metadata['basedir'], '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self.logger.info("Object '{}' loaded with all known data.".format(self.metadata['name']))
        else:
            self.logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            for i in args:
                self.logger.info("Loading {}...".format(i))
                with open(os.path.join(self.metadata['basedir'], '{}/{}.pkl'.format(i, i)), 'rb') as f:
                    self.analysis[i] = cPickle.load(f)
            self.logger.info("Object '{}' loaded with selected data.".format(self.metadata['name']))

    def unload(self, *args):
        """Load data instances into object.

        If 'all' is in argument list, every loaded dataset is unloaded.

        :Arguments:
            *args*
                datasets to unload
        """
        if 'all' in args:
            self.analysis.clear()
            self.logger.info("Object '{}' unloaded of all data.".format(self.metadata['name']))
        else:
            self.logger.info("Unloading selected data from object {}...".format(self.metadata['name']))
            for i in args:
                self.logger.info("Unloading {}...".format(i))
                self.analysis.pop(i, None)
            self.logger.info("Object '{}' unloaded of selected data.".format(self.metadata['name']))

    def save(self):
        """Save base object metadata.

        """
        self._makedirs(self.metadata["basedir"])

        with open(self.metadata['metafile'], 'w') as f:
            yaml.dump(self.metadata, f)

    def _location(self, trajpath, analysisdir):
        """Build analysis directory path from trajectory path.
    
        :Arguments:
            *trajpath*
                path to trajectory
        """
        p = os.path.abspath(trajpath)
        ind = p.split('/').index('WORK')
        seg = '/'.join(p.split('/')[ind+1:])
        seg = os.path.dirname(seg)
    
        return os.path.join(analysisdir, seg)

    def _makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)

    def _load_base(self):
        """Load base object.
        
        """
        with open(self.metadata['metafile'], 'r') as f:
            self.metadata = yaml.load(f)
        self.universe = MDAnalysis.Universe(self.metadata['structure_file'], self.metadata['trajectory_file'])
    
    def _build_metadata(self, **kwargs):
        """Build metadata. Runs on object generation. 
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
        """
        # building core items
        attributes = {'name': kwargs.pop('name', os.path.basename(os.path.dirname(self.metadata['trajectory_file']))),
                      'logfile': os.path.join(self.metadata['basedir'], '{}.log'.format(self.__name__)),
                      'analysis_list': [],
                      'type': self.__name__,
                      }

        for key in attributes.keys():
            if not key in self.metadata:
                self.metadata[key] = attributes[key]

    def _build_attributes(self):
        """Build attributes. Needed each time object is generated.

        """
        # set up logging
        self.logger = logging.getLogger('{}.{}'.format(self.__name__, self.metadata['name']))
        ch = logging.StreamHandler(sys.stdout)
        fh = logging.FileHandler(self.metadata['logfile'])
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)
        self.logger.setLevel(logging.INFO)

class SimSet(object):
    """Base class for a set of simulation objects.

    """
