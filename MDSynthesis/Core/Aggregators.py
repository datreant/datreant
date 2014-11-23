"""
Aggregators are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Sim) in ways
that are user-friendly.

An Aggregator is designed to be user friendly on its own, but it can be used as
a backend by a Container, too.

"""
import Files
import Workers
import MDSynthesis.Containers
from MDAnalysis import Universe
import os
from functools import wraps

class Aggregator(object):
    """Core functionality for information aggregators.

    """
    def __init__(self, container, containerfile, logger):
        """Initialize thin class with references it needs to perform its function.

        """
        self._container = container
        self._containerfile = containerfile
        self._logger = logger

        self._create()

    def _create(self):
        """Initialize object attributes.
        
        """
        pass

class Info(Aggregator):
    """Interface for accessing metadata and status information.

    """
    def __init__(self):
        """
        """
        self.name = None

class SimInfo(Info):
    """Sim-specific bindings.

    """

    #TODO: not updated yet!
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

class GroupInfo(Info):
    """Group-specific bindings.

    """

class Tags(Aggregator):
    """Interface to tags.

    """
    def __repr__(self):
        tags = self.list()

        if not tags:
            out = "No Tags"
        else:
            out = "Tags:"
            for tag in tags:
                out = out + "\t{}\n".format(tag)
        return out

    def __call__(self):
        """Get all tags for the Container as a list.
    
        :Returns:
            *tags*
                list of all tags
        """
        return self._containerfile.get_tags()

    def list(self):
        """Get all tags for the Container as a list.
    
        :Returns:
            *tags*
                list of all tags
        """
        return self._containerfile.get_tags()

    def add(self, *tags):
        """Add any number of tags to the Container.
    
        Tags are individual strings that serve to differentiate Containers from
        one another. Sometimes preferable to categories.
    
        :Arguments:
           *tags*
              Tags to add. Must be convertable to strings using the str() builtin.

        """
        self._containerfile.add_tags(*tags)
    
    def remove(self, *tags, **kwargs):
        """Remove tags from Container.
        
        Any number of tags can be given as arguments, and these will be
        deleted.

        :Arguments:
            *tags*
                Tags to delete.

        :Keywords:
            *all*
                When True, delete all tags [``False``]
        """
        self._containerfile.del_tags(*tags, **kwargs)

class Categories(Aggregator):
    """Interface to categories.

    """
    def __repr__(self):
        categories = self.dict()

        if not categories:
            out = "No Categories"
        else:
            out = "Categories:"
            for key in categories.keys():
                out = out + "\t{}\t: {}\n\t".format(key, categories[key])
        return out

    def __call__(self):
        """Get all categories for the Container as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories 
        """
        return self._containerfile.get_categories()

    def dict(self):
        """Get all categories for the Container as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories 
        """
        return self._containerfile.get_categories()

    def add(self, **categories):
        """Add any number of categories to the Container.

        Categories are key-value pairs of strings that serve to differentiate
        Containers from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will replace
        the value for that category.

        :Keywords:
            *categories*
                Categories to add. Keyword used as key, value used as value. Both
                must be convertible to strings using the str() builtin.
        
        """
        self._containerfile.add_categories(**categories)
    
    def remove(self, *categories, **kwargs):
        """Remove categories from  Container.
        
        Any number of categories (keys) can be given as arguments, and these
        keys (with their values) will be deleted.
         
        :Arguments:
            *categories*
                Categories to delete.

        :Keywords:
            *all*
                When True, delete all categories [``False``]
    
        """
        self._containerfile.del_categories(*categories, **kwargs)

class Universes(Aggregator):
    """Interface to universes.

    """
    def __repr__(self):
        universes = self.list()

        if not tags:
            out = "No Universes"
        else:
            out = "Tags:"
            for tag in tags:
                out = out + "\t{}\n".format(tag)
        return out

    def add(self, handle, topology, *trajectory):
        """Add a universe definition to the Sim object.

        A Universe is an MDAnalysis object that gives access to the details
        of a simulation trajectory. A Sim object can contain multiple universe
        definitions (topology and trajectory pairs), since it is often
        convenient to have different post-processed versions of the same
        raw trajectory.

        :Arguments:
            *handle*
                given name for selecting the universe
            *topology*
                path to the topology file
            *trajectory*
                path to the trajectory file; multiple files may be given
                and these will be used in order as frames for the trajectory

        """
        self._containerfile.add_universe(handle, topology, *trajectory)
    
    def remove(self, handle):
        """Remove a universe definition.
    
        Also removes any selections associated with the universe.

        :Arguments:
            *handle*
                name of universe to delete
        """
        self._containerfile.del_universe(handle)
    
        if self._container._uname == handle:
            self.detach()
    
    #TODO: implement in SimFile
    def rename(self, handle, new_handle):
        """Change the handle used for the specified Universe.
    
        """
        pass

    #TODO: not updated yet!
    def cache(self, handle, location):
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

        self.universe = Universe(structure_c, *trajectory_c)
        self._logger.info("Universe '{}' now cached.".format(self._uname))
    
    #TODO: not updated yet!
    # ensure that trajectories are decached on object destruction
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
        
class Selections(Aggregator):
    """Selection manager for Sims.

    Selections are accessible as items using their handles. Each time they are
    called, they are regenerated from the universe that is currently active. In this
    way, changes in the universe topology are reflected in the selections.

    """
    def __getitem__(self, handle):
        """Get AtomGroup corresponding to the given named selection.
        
        """
        selstring = self._containerfile.get_selection(self._container._uname, handle)
        return self._container.universe.selectAtoms(*selstring)

    def __setitem__(self, handle, selection):
        """Selection for the given handle and the active universe.
    
        """
        self._containerfile.add_selection(self._container._uname, handle, *selection)
    
    def __delitem__(self, handle):
        """Remove stored selection for given handle and the active universe.
    
        """
        self._containerfile.del_selection(self._container._uname, handle)

    def add(self, handle, *selection):
        """Add an atom selection for the attached Universe.
        
        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        :Arguments:
            *handle*
                name to use for the selection
            *selection*
                selection string; multiple strings may be given and their
                order will be preserved, which is useful for e.g. structural 
                alignments
        """
        self._containerfile.add_selection(self._container._uname, handle, *selection)
        
    def remove(self, handle):
        """Remove an atom selection for the attached Universe.
        
        :Arguments:
            *handle*
                name of selection to remove
        """
        self._containerfile.del_selection(self._container._uname, handle)
    
    #TODO: handle case where no universe is loaded
    def keys(self):
        """Return a list of all selection handles.
    
        """
        return self._containerfile.list_selections(self._container._uname)

class Members(Aggregator):
    """Member manager for Groups.

    """
    def __call__(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        """
        return list(self._members)

    def _create(self):
        """Load existing members.

        """
        #TODO: need to route these through Finder.
        self._members = list()
        for uuid in self._containerfile.get_members_uuid():
            member = self._containerfile.get_member(uuid)

            if member['containertype'] == 'Sim':
                self._members.append(MDSynthesis.Containers.Sim(member['abspath'], detached=True))
            elif member['containertype'] == 'Group':
                self._members.append(MDSynthesis.Containers.Group(member['abspath'], detached=True))

    def list():
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        """
        return list(self._members)

    def add(self, *containers):
        """Add any number of members to the Group.

        :Arguments:
            *containers*
                Sims and/or Groups to be added
        """
        for container in containers:
            self._containerfile.add_member(container.info.uuid, container.info.containertype, container.info.location)
    
    def remove(self, *uuids):
        """Remove any number of members from the Group.
    
        :Arguments:
            *uuids*
                the uuids of the members to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        self._containerfile.del_member(*uuids)

class Data(Aggregator):
    """Interface to stored data.

    """
    def __repr__(self):
        data = self.list()

        if not data:
            out = "No Data"
        else:
            out = "Data:"
            for datum in data:
                out = out + "\t{}\n".format(datum)
        return out

    def _makedirs(self, p):
        """Make directories and all parents necessary.

        :Arguments:
            *p*
                directory path to make
        """
        if not os.path.exists(p):
            os.makedirs(p)

    def _get_datafile(self, handle):
        """Return path to datafile corresponding to given handle.

        :Arguments:
            *handle*
                name of dataset whose datafile path to return

        :Returns:
            *datafile*
                datafile path

        """
        datafile = os.path.join(self._containerfile.get_location(), 
                                handle, Files.datafile) 
        return datafile

    #TODO: fundamentally flawed. File is generated if directory exists; revisit
    # immediately
    def _read_datafile(func):
        """Decorator for generating DataFile instance for reading data.

        DataFile instance is generated and mounted at self._datafile. It
        is then dereferenced after the method call. Since data files can be
        deleted in the filesystem, this should handle cleanly the scenarios
        in which data appears, goes missing, etc. while a Container is loaded.

        ``Note``: methods wrapped with this decorator need to have *handle*
        as the first argument.
        
        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            filename = os.path.join(self._containerfile.get_location(), 
                                    handle, Files.datafile) 

            if os.path.exists(filename):
                self._datafile = Files.DataFile(filename, logger=self._logger)
                try:
                    out = func(self, handle, *args, **kwargs)
                finally:
                    del self._datafile
            else:
                self._logger.warning("No data named '{}' present.".format(handle))
                out = None

            return out

        return inner

    def _write_datafile(func):
        """Decorator for generating DataFile instance for writing data.

        DataFile instance is generated and mounted at self._datafile. It
        is then dereferenced after the method call. Since data files can be
        deleted in the filesystem, this should handle cleanly the scenarios
        in which data appears, goes missing, etc. while a Container is loaded.

        ``Note``: methods wrapped with this decorator need to have *handle*
        as the first argument.
        
        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            dirname = os.path.join(self._containerfile.get_location(), handle)
            filename = os.path.join(dirname, Files.datafile) 

            self._makedirs(dirname)
            self._datafile = Files.DataFile(filename, logger=self._logger)

            try:
                out = func(self, handle, *args, **kwargs)
            finally:
                del self._datafile

            return out

        return inner

    def __getitem__(self, handle):
        """Get dataset corresponding to given handle.

        If dataset doesn't exist, ``None`` is returned.
        
        :Arguments:
            *handle*
                name of data to retrieve

        :Returns:
            *data*
                stored data; ``None`` if nonexistent
        """
        return self.retrieve(handle)

    def __setitem__(self, handle, data):
        """Set dataset corresponding to given handle.
        
        A data instance must be either a pandas Series, DataFrame, or Panel
        object. If dataset doesn't exist, it is added. If a dataset already
        exists for the given handle, it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data to store; must be a pandas Series, DataFrame, or Panel

        """
        self.add(handle, data)
    
    def __delitem__(self, handle):
        """Remove a dataset.

        Note: the directory containing the dataset file (``Data.h5``) will NOT
        be removed if it still contains file after the removal of the dataset
        file.

        :Arguments:
            *handle*
                name of dataset to delete
    
        """
        self.remove(handle)

    @_write_datafile
    def add(self, handle, data):
        """Store data in Container.

        A data instance must be either a pandas Series, DataFrame, or Panel
        object. If dataset doesn't exist, it is added. If a dataset already
        exists for the given handle, it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data to store; must be a pandas Series, DataFrame, or Panel

        """
        self._datafile.add_data('main', data)
        
    def remove(self, handle, **kwargs):
        """Remove a dataset, or some subset of a dataset.

        Note: in the case the whole dataset is removed, the directory
        containing the dataset file (``Data.h5``) will NOT be removed if it
        still contains file(s) after the removal of the dataset file.

        :Arguments:
            *handle*
                name of dataset to delete
            *where*
                conditions for what rows/columns to remove
            *start* 
                row number to start selection
            *stop*  
                row number to stop selection
    
        """

        if kwargs: 
            if os.path.exists(self._get_datafile(handle)):
                self._delete_data(handle, **kwargs)
            else:
                self._logger.info("No data named '{}' present".format(handle))
        else:
            datafile = self._get_datafile(handle)

            if os.path.exists(self._get_datafile(handle)):
                os.remove(datafile)
                try: 
                    os.rmdir(os.path.dirname(datafile))
                except OSError:
                    self._logger.info("Directory for '{}' not empty. Retaining directory.".format(handle))
            else:
                self._logger.info("No data named '{}' present. Nothing to remove".format(handle))

    @_write_datafile
    def _delete_data(self, handle, **kwargs):
        """Remove a dataset, or some subset of a dataset.

        This method loads the given data instance before doing anything,
        and should generally not be used to remove data since it will create
        a datafile object if one is not already present, which could have
        side-effects for other instances of this Container.

        Note: in the case the whole dataset is removed, the directory
        containing the dataset file (``Data.h5``) will NOT be removed if it
        still contains file(s) after the removal of the dataset file.

        :Arguments:
            *handle*
                name of dataset to delete
            *where*
                conditions for what rows/columns to remove
            *start* 
                row number to start selection
            *stop*  
                row number to stop selection

        """
        try:
            self._datafile.del_data('main', **kwargs)
        except NotImplementedError:
            self._logger.info("Dataset '{}' empty; removing.".format(handle))
            datafile = self._get_datafile(handle)

            if os.path.exists(self._get_datafile(handle)):
                os.remove(datafile)
                try: 
                    os.rmdir(os.path.dirname(datafile))
                except OSError:
                    self._logger.info("Directory for '{}' not empty. Retaining directory.".format(handle))
            else:
                self._logger.info("No data named '{}' present. Nothing to remove".format(handle))

    @_read_datafile
    def retrieve(self, handle, **kwargs):
        """Retrieve stored data.

        The pandas object (Series, DataFrame, or Panel) for the given handle
        is returned.
    
        If dataset doesn't exist, ``None`` is returned.

        Subsets of the whole dataset can be returned using keywords such as
        *start* and *stop* for ranges of rows, and *columns* for selected
        columns.

        The *where* keyword takes a string as input and can be used to filter
        out rows and columns without loading the full object into memory. For
        example, given a DataFrame with handle 'mydata' with columns (A, B, C,
        D), one could return all rows for columns A and C for which column D is
        greater than .3 with::

            retrieve('mydata', where='columns=[A,C] & D > .3')

        Or, if I wanted all rows with index = 3 (there could be more than
        one)::

            retrieve('mydata', where='index = 3')

        See :meth:pandas.HDFStore.select() for more information.
        
        :Arguments:
            *handle*
                name of data to retrieve

        :Keywords:
            *where*
                conditions for what rows/columns to return
            *start* 
                row number to start selection
            *stop*  
                row number to stop selection
            *columns*
                list of columns to return; all columns returned by default
            *iterator* 
                if True, return an iterator [``False``]
            *chunksize*
                number of rows to include in iteration; implies
                ``iterator=True`` 

        :Returns:
            *data*
                stored data; ``None`` if nonexistent

        """
        return self._datafile.get_data('main', **kwargs)

    @_write_datafile
    def append(self, handle, data):
        """Append rows to an existing dataset.
    
        The object must be of the same class (Series, DataFrame, Panel) as the
        existing dataset, and it must have exactly the same columns (names
        included).

        :Arguments:
            *handle*
                name of data to append to
            *data*
                data to append

        """
        self._datafile.append_data('main', data)
    
    def list(self):
        """List available datasets.

        :Returns:
            *handles*
                list of handles to available datasets

        """
        datasets = list()
        for root, dirs, files in os.walk(self._containerfile.get_location()):
            if Files.datafile in files:
                datasets.append(os.path.basename(root))

        datasets.sort()
    
        return datasets

class Bunch(object):
    def __init__(self, odict):
        adict = dict(odict)
        for key in adict:
            if type(adict[key]) is dict:
                adict[key] = RwBunch(adict[key])
        self.__dict__ = adict

class Database(Aggregator):
    """Database object for tracking and coordinating Containers.

    The Database object stores information on all Containers it is made aware of.
    This centralized storage allows Containers to find each other when necessary;
    this is especially important for Groups.

    This object is the interface of Container objects to the database file.

    """

    def __init__(self, database, **kwargs):
        """Generate Database object for the first time, or interface with an existing one.

        :Arguments:
            *database*
                directory containing a Database file; if no Database file is
                found, a new one is generated

        """
        super(Database, self).__init__()
        self._database = dict()              # the database data itself

        database = os.path.abspath(database)
        dbfile = os.path.join(database, self._databasefile)
        if os.path.exists(dbfile):
            self._regenerate(database, **kwargs)
        else:
            self._generate(database, **kwargs)

    def _generate(self, database):
        """Generate a new database.

        """
        self._database['basedir'] = database
        self._build_metadata()
        self._build_attributes()

        # write to database file
        self.commit()
        self._start_logger(database)

    def _regenerate(self, database):
        """Re-generate existing database.

        """
        self.database['basedir'] = database
        self.refresh()
        self._start_logger(database)

        self._check_location(database)

        # rebuild missing parts
        self._build_metadata()
        self._build_attributes()

    def _handshake(self):
        """Run check to ensure that database is fine.

        """
        #TODO: add various checks to ensure things are in working order
        return ('basedir' in self.database)

    def search(self, searchstring):
        """Search the Database for Containers that match certain criteria.

        Results are printed in digest form to the screen. To print full
        metadata for all matching containers, use print='full'

        :Arguments:
            *searchstring*
                string giving the search query

        :Keywords:
            *print*
                format of results printed to ouptut

        :Returns:
            *locations*
                filesystem paths to Containers that match criteria

        """
        #TODO: Add in selection system similar to that implemented in
        # MDAnalysis for atom selections. This one, however, will parse
        # metadata elements, and shouldn't be quite so complex

        return

    def add(self, *containers, **kwargs):
        """Add Container to Database.

        :Arguments:
            *containers*
                Containers to add, each given as a path to a Container directory
                or as a generated Container object

        """
        for container in containers:
            if isinstance(container, basestring) and os.path.isdir(container):
                with self.util.open(os.path.join(container, self._containerfile), 'r') as f:
                    meta = yaml.load(f)
                uuid = meta['uuid']
                meta['basedir'] = os.path.abspath(container)
                self.database['containers'][uuid] = meta
                with self.util.open(os.path.join(container, self._containerfile), 'w') as f:
                    yaml.dump(meta, f)
            else:
                uuid = container.metadata['uuid']
                self.database['containers'][uuid] = container.metadata

            self.database['containers'][uuid]['database'] = self.database['basedir']

            # since this method is used for Container init, basedir may not
            # be defined in metadata yet
            if not ('basedir' in self.database['containers'][uuid]):
                container.metadata['basedir'] = self._build_basedir(uuid)
                self.database['containers'][uuid]['basedir'] = container.metadata['basedir']
                with self.util.open(os.path.join(container.metadata['basedir'], self._containerfile), 'w') as f:
                    yaml.dump(self.database['containers'][uuid], f)
                self.commit()
            else:
                self.push(uuid)
            self._logger.info("Added {} container '{}' to database.".format(self.database['containers'][uuid]['class'], self.database['containers'][uuid]['name']))

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
            *all*
                if True, will remove all entries ``[False]``
        """
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['containers'] ]

        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            matches = []
            for entry in self.database['containers'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])

            for match in matches:
                self.database['containers'].pop(match, None)

    def clean(self):
        """Clear entries from Database corresponding to Containers that can't be found.

        """
        self._logger.info("Cleaning out entries that cannot be found.")

        uuids = [ x for x in self.database['containers'] ]
        self._get_containers(*uuids)

        for uuid in uuids:
            if not self.database['containers'][uuid]['basedir']:
                self._logger.info("Removing: {} ({})".format(self.database['containers'][uuid]['name'], uuid))
                self.database['containers'].pop(uuid)
        self._logger.info("Database is clean.")

    def commit(self):
        """Save the current state of the database to its file.

        """
        self.util.makedirs(self.database['basedir'])
        with self.util.open(os.path.join(self.database['basedir'], self._databasefile), 'w') as f:
            yaml.dump(self.database, f)

    def refresh(self):
        """Reload contents of database file.

        """
        dbfile = os.path.join(self.database['basedir'], self._databasefile)
        with self.util.open(dbfile, 'r') as f:
            self.database = yaml.load(f)

    def pull(self, *containers, **kwargs):
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
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['containers'] ]

        matches = []
        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            for entry in self.database['containers'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])

        # ensure we are finding the right Container
        basedirs = self._get_containers(*matches)

        for i in xrange(len(matches)):
            if basedirs[i]:
                with self.util.open(os.path.join(basedirs[i], self._containerfile), 'r') as f:
                    self.database['containers'][matches[i]] = yaml.load(f)
        self.commit()

    def push(self, *containers, **kwargs):
        """Update Container metadata with information stored in Database.

        This is the opposite of `:meth:self.pull()`

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
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['containers'] ]

        matches = []
        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            for entry in self.database['containers'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])

        # since this method is used for Container init, basedir may not
        # be defined in metadata yet
        for match in matches:
            if not ('basedir' in self.database['containers'][match]):
                self.database['containers'][match]['basedir'] = self._build_basedir(match)

        # ensure we are finding the right Container
        basedirs = self._get_containers(*matches)

        for i in xrange(len(matches)):
            if basedirs[i]:
                with self.util.open(os.path.join(basedirs[i], self._containerfile), 'w') as f:
                    yaml.dump(self.database['containers'][matches[i]], f)
        self.commit()

    def _get_containers(self, *uuids):
        """Get path to Containers.

        Will perform checks to ensure the Container returned matches the uuid given.
        It will go looking for the Container if not found at last known location.

        :Arguments:
            *uuids*
                unique ids for Containers to return

        :Returns:
            *containers*
                tuple giving paths to Containers
        """
        containers = [None]*len(uuids)
        missing = [None]*len(uuids)
        for i in xrange(len(uuids)):
            if not self.database['containers'][uuids[i]]['basedir']:
                continue

            if os.path.exists(os.path.join(self.database['containers'][uuids[i]]['basedir'], self._containerfile)):
                with self.util.open(os.path.join(self.database['containers'][uuids[i]]['basedir'], self._containerfile), 'r') as f:
                    meta = yaml.load(f)
                if meta['uuid'] == uuids[i]:
                    containers[i] = self.database['containers'][uuids[i]]['basedir']
                else:
                    self._logger.info("Missing: {} ({})".format(self.database['containers'][uuids[i]]['name'], uuids[i]))
                    missing[i] = uuids[i]
            else:
                self._logger.info("Missing: {} ({})".format(self.database['containers'][uuids[i]]['name'], uuids[i]))
                missing[i] = uuids[i]

        if any(missing):
            missing = self._locate_containers(*missing)

        # build final list of paths
        for i in xrange(len(uuids)):
            if not containers[i]:
                containers[i] = missing[i]

        return containers

    def discover(self):
        """Traverse filesystem downward from Database directory and add all new Containers found.

        """
        for root, dirs, files in os.walk(self.database['basedir']):
            if self._containerfile in files:
                dirs = []
                self.add(root)
        self.commit()

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

    def _check_location(self, database, **kwargs):
        """Check Database location; if changed, send new location to all Containers.

        :Keywords:
            *force*
                if True, new location sent to all Containers even if unchanged;
                default False
        """
        force = kwargs.pop('force', False)
        database = os.path.abspath(database)

        if (database != self.database['basedir']) or force:
            self.database['basedir'] = database

            # update entries first
            self.pull(all=True)

            for entry in self.database['containers'].values():
                entry['database'] = self.database['basedir']

            self.commit()
            self.push(all=True)

    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.

        Only adds keys; never modifies existing ones.

        """
        attributes = {'class': self.__class__.__name__,
                      'name': kwargs.pop('name', os.path.basename(self.database['basedir'])),
                      'containers': dict(),
                      }

        for key in attributes:
            if not key in self.database:
                self.database[key] = attributes[key]

    def _build_attributes(self):
        """Build attributes of Database. Runs each time object is generated.

        """

    def _locate_containers(self, *uuids):
        """Find Containers by traversing downward through the filesystem.

        Looks in each directory below the Database. If found, the basedir for the
        Container is updated in both metadata and the Database.

        :Arguments:
            *uuids*
                unique ids for Containers to return
        """
        self._logger.info("Searching for {} Containers.".format(len(uuids) - uuids.count(None)))
        containers = [None]*len(uuids)
        for root, dirs, files in os.walk(self.database['basedir']):
            if self._containerfile in files:
                dirs = []
                with self.util.open(os.path.join(root, self._containerfile), 'r') as f:
                    meta = yaml.load(f)
                try:
                    i = uuids.index(meta['uuid'])
                    containers[i] = os.path.abspath(root)
                    meta['basedir'] = containers[i]

                    # update basedir in Container metadata and in Database
                    with self.util.open(os.path.join(root, self._containerfile), 'w') as f:
                        yaml.dump(meta, f)
                    self.database['containers'][uuids[i]]['basedir'] = containers[i]
                    self._logger.info("Found: {} ({})\nLocation: {}".format(meta['name'], uuids[i], meta['basedir']))
                except ValueError:
                    pass

        for i in xrange(len(containers)):
            if uuids[i]:
                if not containers[i]:
                    self.database['containers'][uuids[i]]['basedir'] = None
                    self._logger.warning("Not found: {} ({})".format(self.database['containers'][uuids[i]]['name'], uuids[i]))

        self._logger.info("{} Containers not found.".format(containers.count(None) - uuids.count(None)))
        return containers

    def _build_basedir(self, uuid):
        """Build basedir location based on database location, Container class, and Container name.

        :Arguments:
            *database*
                directory where database resides
            *name*
        """
        database = self.database['basedir']
        container = self.database['containers'][uuid]

        # if name given and directory with name doesn't already exist, make named basedir
        if container['name'] and not os.path.exists(os.path.join(database, container['class'], container['name'])):
            dest = container['name']
        # if basedir already exists, use UUID instead
        else:
            dest = container['uuid']

        dest = os.path.join(database, container['class'], dest)
        self.util.makedirs(dest)

        return dest

    def _start_logger(self, basedir):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.database['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)

            # file handler
            logfile = os.path.join(basedir, self._containerlog)
            fh = logging.FileHandler(logfile)
            ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
            fh.setFormatter(ff)
            self._logger.addHandler(fh)

            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)

