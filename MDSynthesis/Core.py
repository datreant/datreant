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
    
    def save(self, *args):
        """Save base object metadata and, if desired, analysis data instances.

        By providing names of loaded datasets as arguments, you can save the
        loaded versions of the data to their source files. This is useful if
        you have need to manually edit the data in a source file, as you can
        load it, make changes, then write it out.

        If no arguments are given, then no datasets are saved to their source
        files. Only metadata is saved.

        If 'all' is in argument list, every dataset that is loaded is written
        to its file.

        :Arguments:
            *args*
                datasets to load
            
        :Keywords:
            *force*
                if True, reload data even if already loaded; default False

        """
        basedir = self._rel2abspath(self.metadata["basedir"])
        self._makedirs(basedir)

        with open(os.path.join(basedir, self.metadata['metafile']), 'w') as f:
            yaml.dump(self.metadata, f)

        if 'all' in args:
            self._logger.info("Saving all loaded data into source files for '{}'...".format(self.metadata['name']))
            savelist = self.analysis
        elif len(args) != 0:
            self._logger.info("Saving selected data into source files for '{}'...".format(self.metadata['name']))
            savelist = args

        for i in savelist:
            self._logger.info("Saving {}...".format(i))
            with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}.pkl'.format(i, i)), 'wb') as f:
                cPickle.dump(self.analysis[i], f)
        self._logger.info("All selected data saved.")

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
            loadlist = self.metadata['analyses']
        else:
            self._logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            loadlist = args

        for i in loadlist:
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

    def _update_database(self):
        """Update metadata stored in Database for this Container.
    
        """

    def _update_container(self):
        """Update metadata from information stored in Database for this Container.

        Note: This will overwrite metadata file with Database version!

        """

    def _locate_database(self):
        """Find a database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. If it does not find
        it, it creates a new one where it things it should have been.

        """
    
    def _generate_database(self):
        """Generate a database for the first time.

        If no Database exists for the Container to give its information to,
        this method is used to create a new one.

        """

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

class Database(object):
    """Database object for tracking and coordinating Containers.
    
    """
    def __init__(self):
        """Generate Database object for the first time, or interface with an existing one.


        """
    
    def search(self):
        """Search the Database for Containers that match certain criteria.

        Results are printed in digest form to the screen. To print full
        metadata for all matching containers, use print='full'

        :Keywords:
            *print*
                format of results printed to ouptut

        :Returns:
            *locations*
                filesystem paths to Containers that match criteria

        """

    def add(self, *containers, **kwargs):
        """Add Container to Database.

        :Arguments:
            *containers*
                Containers to add, each given as a path to a Container directory
                or as a generated Container object
            
        """

    def remove(self, *containers, **kwargs):
        """Remove Container from Database.

        Note: if Container name is used to specify removal and more than one
        Container has that name, then both will be removed.

        :Arguments:
            *containers*
                Containers to remove, each given as a path to a Container directory,
                a Container UUID, or a Container's given name

        :Keywords:
            *hard*
                if True, delete Container object from filesystem too ``[False]``

        """

    def clean(self):
        """Clear entries from Database corresponding to Containers that can't be found.

        """

    def update_database(self, *args, **kwargs):
        """Update information stored in Database from Container metadata.

        Note: if Container name is used to specify the update, all Containers
        with that name will be updated in the Database.

        :Arguments:
            *args*
                Containers to update, each given as a path to a Container directory,
                a Container UUID, or a Container's given name

        :Keywords:
            *all*
                if True, will update entries for all known Containers from metadata files
        """

    def update_container(self, *containers):
        """Update Container metadata with information stored in Database.

        This is the opposite of `:meth:self.update_database()`

        Note: if Container name is used to specify the update, all Containers
        with that name will have metadata updated.

        :Arguments:
            *containers*
                Containers to update; either a path to a Container directory,
                Container UUID, or a Container's given name
        :Keywords:
            *all*
                if True, will update all known Container metadata files from entries

        """

    def discover_containers(self):
        """Traverse filesystem downward from Database directory and add all new Containers found.
        
        """
    
    def _check_location(self):
        """Check Database location; if changed, send new location to all Containers.

        """
        # update projectdir

    def _locate_container(self):
        """Find a Container that has moved by traversing downward through the filesystem. 
            
        """

