"""
Lower level mixins.

"""
import os
import yaml
import sys
import cPickle
import logging
from multiprocessing import Process

class ContainerCore(object):
    """Mixin class for all Containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """
    def __init__(self):
        """
        
        """
        self.metadata = dict()              # information about object; defines base object
        self.analysis = dict()              # analysis data 'modular dock'
    
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

    def load(self, *args, **kwargs):
        """Load data instances into object.

        If 'all' is in argument list, every available dataset is loaded.

        :Arguments:
            *args*
                datasets to load
            
        :Keywords:
            *force*
                if True, reload data even if already loaded; default False
        """

        force = kwargs.pop('force', False)

        if 'all' in args:
            self._logger.info("Loading all known data into object '{}'...".format(self.metadata['name']))
            for i in self.metadata['analysis_list']:
                if (i not in self.analysis) or (force == True):
                    self._logger.info("Loading {}...".format(i))
                    with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                        self.analysis[i] = cPickle.load(f)
                else:
                    self._logger.info("Skipping reload of {}...".format(i))
            self._logger.info("Object '{}' loaded with all known data.".format(self.metadata['name']))
        else:
            self._logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            for i in args:
                if (i not in self.analysis) or (force == True):
                    self._logger.info("Loading {}...".format(i))
                    with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'rb') as f:
                        self.analysis[i] = cPickle.load(f)
                else:
                    self._logger.info("Skipping reload of {}...".format(i))
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

    def _abs2relpath(self, abspath):
        """Return path to file relative to project path.
        
        """
        return abspath.replace(self.metadata['projectdir'], '$PROJECT')

    def _rel2abspath(self, relpath):
        """Return realpath given a path relative to project directory. The
            opposite of _project_relpath.
        """
        return relpath.replace('$PROJECT', self.metadata['projectdir'])
    
    def _update_projectdir(self, basedir_abs):
        """Update projectdir based on given (and current) absolute path of basedir.

        """
        self.metadata['projectdir'] = basedir_abs.partition('/MDSynthesis')[0]
        
    def _makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default is trajectory file directory
                basename
        """
        # building core items
        attributes = {'name': kwargs.pop('name', self.__class__.__name__),
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

class OperatorCore(object):
    """Mixin class for all Operators.

    The OperatorCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Operator objects.

    """
    def __init__(self, *args, **kwargs):
        """
        
        """
        self.systems = list(args)

    def run(self, **kwargs):
        """Obtain compute-intensive data, usually timeseries.

        :Keywords:
            *force*
                If True, force recollection of data; default False

            **kwargs passed to `:meth:self._run_system()`
        """
        joblist = []
        force = kwargs.pop('force', False)

        for system in self.systems:
            if (not self._datacheck(system)) or force:
                p = (Process(target=self._run_system, args=(system,), kwargs=kwargs))
                p.start()
                joblist.append(p)
            else:
                system._logger.info('{} data already present; skipping data collection.'.format(self.__class__.__name__))

            # update analysis list in each object
            if not self.__class__.__name__ in system.metadata['analysis_list']:
                system.metadata['analysis_list'].append(self.__class__.__name__)
                system.save()

        for p in joblist:
            p.join()

    def analyze(self, **kwargs):
        """Perform analysis of compute-intensive.

        Does not require stepping through any trajectories.

        """
        # make sure data loaded into each system; should use try/catch here
        self._load()

    def _save(self, system, sys_results):
        """Save results to main data file.

        :Arguments:
            *system*
                system to save data for
            *sys_results*
                results for system
        """
        analysis_dir = self._make_savedir(system)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__class__.__name__))

        with open(main_file, 'wb') as f:
            cPickle.dump(sys_results, f)

    def _make_savedir(self, system):
        """Make directory where all output files are placed.

        :Arguments:
            *system*
                system to save data for

        :Returns:
            *analysis_dir*
                full path to output file directory

        """
        analysis_dir = os.path.join(system._rel2abspath(system.metadata['basedir']), self.__class__.__name__)
        system._makedirs(analysis_dir)

        return analysis_dir
    
    def _load(self, **kwargs):
        """Load data for each system if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for system in self.systems:
            if (not self.__class__.__name__ in system.analysis.keys()) or force:
                system.load(self.__class__.__name__)
    
    def _datacheck(self, system):
        """Check if data file already present.

        :Arguments:
            *system*
                Container object
                
        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        analysis_dir = os.path.join(system._rel2abspath(system.metadata['basedir']), self.__class__.__name__)
        main_file = os.path.join(analysis_dir, '{}.pkl'.format(self.__class__.__name__))
        return os.path.isfile(main_file)
