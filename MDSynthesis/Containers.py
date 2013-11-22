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
    """Base class for simulation objects.

    """

    def __init__(self, *args, **kwargs):
        """Generate a Sim object.

        To regenerate an existing Sim object, give as *system* a directory that
        contains a Sim object metadata file (self.__class__.__name__ + ".yaml").

        :Arguments:

        :Keywords:
            NOTE: keywords only used if *system* is a universe.
            *name*
                desired name for object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
            *location*
                where to store object, only used if *system* is a universe;
                default automatically places object in MDSynthesis directory
                structure. See the :mod:MDSynthesis documentation for details.
            *projectdir*
                path to main project directory; required if no *location* given
            *pluck_segment*
                tuple with components of *trajpath* to leave out of final Sim
                object directory path, e.g. ('WORK/',)
        """
        self.metadata = dict()              # information about object; defines base object
        self.selections = dict()            # AtomGroups
        self.analysis = dict()              # analysis data 'modular dock'

        if (os.path.isdir(args[0])):
        # if system is a directory string, load existing object
            self._regenerate(*args)
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
        """Unload data instances from object.

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
        basedir = self._rel2abspath(self.metadata["basedir"])
        self._makedirs(basedir)

        with open(os.path.join(basedir, self.metadata['metafile']), 'w') as f:
            yaml.dump(self.metadata, f)

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
        location = kwargs.pop('location', None)
        if location == None:
            try:
                projectdir = kwargs.pop('projectdir')
            except KeyError:
                raise KeyError("Cannot construct {} object without location or projectdir. See object documentation for details.".format(self.__class__.__name__))
            self.metadata['projectdir'] = os.path.abspath(projectdir)

            # process plucked segments
            pluck_segment = kwargs.pop('pluck_segment', ('',))
            if isinstance(pluck_segment, basestring):
                pluck_segment = (pluck_segment,)

            self.metadata["basedir"] = self._build_location(system.trajectory.filename, *pluck_segment)

        else:
            self.metadata['projectdir'] = os.path.abspath(location)
            self.metadata['basedir'] = '$PROJECT/MDSynthesis/{}/{}'.format(self.__class__.__name__, kwargs.get('name', os.path.splitext(os.path.basename(system.filename))[0]))

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
        # building core items
        attributes = {'name': kwargs.pop('name', os.path.basename(os.path.dirname(self.metadata['trajectory_files'][0]))),
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
        self.logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))
        self.logger.setLevel(logging.INFO)

        # file handler
        logfile = self._rel2abspath(os.path.join(self.metadata['basedir'], self.metadata['logfile']))
        fh = logging.FileHandler(logfile)
        ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        fh.setFormatter(ff)
        self.logger.addHandler(fh)

        # output handler
        ch = logging.StreamHandler(sys.stdout)
        cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        ch.setFormatter(cf)
        self.logger.addHandler(ch)
        
class SimSet(object):
    """Base class for a set of simulation objects.

    """
