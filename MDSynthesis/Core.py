"""
Lower level mixins.

"""
import os
import yaml
import sys
import cPickle
import logging
import glob
from uuid import uuid4
from multiprocessing import Process

metafile = 'metadata.yaml'
logfile = 'logfile.log'
datafile = 'data.pkl'
database = 'MDSdatabase.yaml'

class ContainerCore(object):
    """Mixin class for all Containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """
    self._metafile = metafile
    self._logfile = logfile
    self._datafile = datafile
    self._database = database

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
                datasets to save

        """
        self._makedirs(self.metadata['basedir'])

        with open(os.path.join(basedir, self._metafile), 'w') as f:
            yaml.dump(self.metadata, f)

        if 'all' in args:
            self._logger.info("Saving all loaded data into source files for '{}'...".format(self.metadata['name']))
            savelist = self.analysis
        elif len(args) != 0:
            self._logger.info("Saving selected data into source files for '{}'...".format(self.metadata['name']))
            savelist = args

        for i in savelist:
            self._logger.info("Saving {}...".format(i))
            with open(os.path.join(self._rel2abspath(self.metadata['basedir']), '{}/{}'.format(i, self._datafile)), 'wb') as f:
                cPickle.dump(self.analysis[i], f)
        self._logger.info("All selected data saved.")

    def refresh(self):
        """Reloads metadata from file.

        """
        metafile = os.path.join(self.metadata['basedir'], self._metafile)
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
                with open(os.path.join(self.metadata['basedir'], '{}/{}'.format(i, self._datafile)), 'rb') as f:
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
                object in some analyses; default None
        """
        # building core items
        uuid = self._generate_id()
        attributes = {'id': uuid,
                      'name': kwargs.pop('name', None),
                      'analysis': list(),
                      'class': self.__class__.__name__,
                      'categories': kwargs.pop('categories', dict()),
                      'tags': kwargs.pop('tags', list()),
                      }

        for key in attributes:
            if not key in self.metadata:
                self.metadata[key] = attributes[key]
    
    def _build_basedir(self, database, name):
        """Build basedir location based on database location, Container class, and Container name.

        :Arguments:
            *database*
                directory where database resides
            *name*
        """
        basedir = os.path.join(os.path.abspath(database), '{}/{}'.format(self.__class__.__name__, name)

        if os.path.exists(basedir):
            basedir = self._build_basedir(database, self.metadata['id'])

        




        return os.path.join(os.path.abspath(database), '{}/{}'.format(self.__class__.__name__, name)

    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)

            # file handler
            logfile = os.path.join(self.metadata['basedir'], self._logfile)
            fh = logging.FileHandler(logfile)
            ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
            fh.setFormatter(ff)
            self._logger.addHandler(fh)

            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)

    def _generate_id(self):
        """Generate a 'unique' identifier.

        """
        return str(uuid4())

    def _update_database(self):
        """Update metadata stored in Database for this Container.
    
        """

    def _update_container(self):
        """Update metadata from information stored in Database for this Container.

        Note: This will overwrite metadata file with Database version!

        """


    def _init_database(self, **kwargs):
        """On generation of Container, perform standard database interactions.

        :Keywords:
            *database*
                path to database; default None
            *locate*
                if True, automatically try to find a database if none specified;
                if False, generate new database in current directory if none 
                specified; default True
        
        """
        database = kwargs.pop('database', None)
        locate = kwargs.pop('locate', True)

        if database:
            database = os.path.dirname(database)
            self._logger.info("Attempting to connect to database in {}".format(database))
            self._connect_database(database)
        else:
            self._logger.info("No database specified. Looking upward to find one.".format(database))
            database = self._locate_database(startdir='.')
            if not database:
                self.logger.info("Generating new database in current directory.")
                self._connect_database('.')
            
    def _connect_database(self, **kwargs):
        """Connect Container to a Database.

        If the Database doesn't exist, it will be created.

        :Keywords:
            *database*
                path to database

        :Returns:
            *success*
                True if connection succeeded; False otherwise
        """
        # attempt to open database
        database = kwargs.pop('database', self.metadata['database'])
        db = Database(database)
        
        if db._handshake():
            self._logger.info("Handshake success; database now in {}".format(db.database['basedir'])
            self.metadata['database'] = db.metadata['basedir']
            db.add(self)
            success = True
        else:
            self._logger.warning("Specified database failed handshake; not a real database?")
            success = False
    
        return success

    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The first such file
        found will become its new Database. If none is found, 

        :Keywords:
            *startdir*
                directory from which to begin upward search; default is
                Container basedir

        :Returns:
            *database*
                path to located database; if no database found, is None
        
        """
        startdir = kwargs.pop('startdir', self.metadata['basedir'])

        # search upward for a database
        directory = startdir
        found = False
        
        self._logger.info("Beginning search for database from {}".format(directory))
        while (directory != '/') and (not found):
            directory, tail = os.path.split(directory)
            candidates = glob.glob(os.path.join(directory, self._database))
            
            if candidates:
                self._logger.info("Database candidate located: {}".format(candidates[0])
                found = self._connect_database(candidates[0]):
        
        if not found:
            self._logger.info("No database found!")
            self.metadata['database'] = None

        return self.metadata['database']

class OperatorCore(object):
    """Mixin class for all Operators.

    The OperatorCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Operator objects.

    """
    self._datafile = datafile

    def __init__(self, *args, **kwargs):
        """
        
        """
        self.containers = list(args)

    def run(self, **kwargs):
        """Obtain compute-intensive data, usually timeseries.

        :Keywords:
            *force*
                If True, force recollection of data; default False

            **kwargs passed to `:meth:self._run_container()`
        """
        joblist = []
        force = kwargs.pop('force', False)

        # run analysis on each container as a separate process
        for container in self.containers:
            if (not self._datacheck(system)) or force:
                p = (Process(target=self._run_container, args=(container,), kwargs=kwargs))
                p.start()
                joblist.append(p)
            else:
                system._logger.info('{} data already present; skipping data collection.'.format(self.__class__.__name__))

        for p in joblist:
            p.join()
    
        # finish up
        for container in self.containers:
            # update analysis list in each object
            if not self.__class__.__name__ in container.metadata['analysis_list']:
                container.metadata['analysis'].append(self.__class__.__name__)
                container.save()

    def analyze(self, **kwargs):
        """Perform analysis of compute-intensive data.

        Does not require stepping through any trajectories.

        """
        # make sure data loaded into each container; should use try/catch here
        self._load()

    def _save(self, container, cont_results):
        """Save results to main data file.

        :Arguments:
            *container*
                container to save data for
            *cont_results*
                results for container
        """
        outputdir = self._make_savedir(container)
        main_file = self._datafile(container)

        with open(main_file, 'wb') as f:
            cPickle.dump(cont_results, f)

    def _make_savedir(self, container):
        """Make directory where all output files are placed.

        :Arguments:
            *container*
                Container for which to save data 

        :Returns:
            *outputdir*
                full path to output directory

        """
        outputdir = self._outputdir(container)
        container._makedirs(outputdir)

        return outputdir
    
    def _load(self, **kwargs):
        """Load data for each container if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for container in self.containers:
            if (not self.__class__.__name__ in container.analysis.keys()) or force:
                container.load(self.__class__.__name__)
    
    def _datacheck(self, container):
        """Check if data file already present.

        :Arguments:
            *container*
                Container object to check
                
        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        outputdir = self._outputdir(container)
        main_file = os.path.join(analysisdir, self._datafile)
        return os.path.isfile(main_file)

    def _outputdir(self, container):
        """Return path to output directory for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *analysis_path*
                path to output directory
        """
        return os.path.join(container.metadata['basedir'], self.__class__.__name__)
    
    def _datafile(self, container):
        """Return path to main output datafile for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *datafile_path*
                path to datafile
        """
        return os.path.join(self._outputdir(container), self._datafile)

class Database(object):
    """Database object for tracking and coordinating Containers.
    
    """
    #TODO: add simple lock scheme to Database so that only one instance at a
    # time can commit changes to the file, and that when this happens all
    # Database instances have a way of knowing a change has been made and can
    # update their record 

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

    def select(self, *containers):

    def deselect(self, *containers):

    def add(self, *containers, **kwargs):
        """Add Container to Database.

        :Arguments:
            *containers*
                Containers to add, each given as a path to a Container directory
                or as a generated Container object
            
        """
        for container in containers:
            self.

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
        # use os.walk
    
    def merge(self, database):
        """Merge another database's contents into this one.
        
        :Arguments:
            *database*
                path to database or Database object

        """

    def split(self, database):
        """Split selected Containers off of database into another.
        
        :Arguments:
            *database*
                path to destination database or Database object
        """
    def _check_location(self):
        """Check Database location; if changed, send new location to all Containers.

        """
        # update projectdir

    def _locate_container(self):
        """Find a Container that has moved by traversing downward through the filesystem. 
            
        """
        # use os.walk

    def _handshake(self):
        """Run check to ensure that database is fine. Return database ID.

        """
        #TODO: add various checks to ensure things are in working order
        return self.database['metadata']['id']

