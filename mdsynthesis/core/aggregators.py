"""
Aggregators are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Sim) in ways
that are user-friendly.

In short, an Aggregator is designed to be user friendly on its own, but it can
be used as a backend by a Container, too.

"""
from MDAnalysis import Universe
import os
from functools import wraps

import persistence
import filesystem
import mdsynthesis as mds

class Aggregator(object):
    """Core functionality for information aggregators.

    """
    def __init__(self, container, containerfile, logger):
        self._container = container
        self._containerfile = containerfile
        self._logger = logger

class Tags(Aggregator):
    """Interface to tags.

    """
    def __repr__(self):
        return "<Tags({})>".format(self.list())

    def __str__(self):
        tags = self.list()
        agg = "Tags"
        majsep = "="
        seplength = len(agg)

        if not tags:
            out = "No Tags"
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for i in xrange(len(tags)):
                out = out + "'{}'\n".format(tags[i])
        return out

    def __iter__(self):
        return self._containerfile.get_tags().__iter__()

    def __len__(self):
        return len(self._containerfile.get_tags())

    def list(self):
        """Get all tags for the Container as a list.
    
        :Returns:
            *tags*
                list of all tags
        """
        tags = self._containerfile.get_tags()
        tags.sort()
        return tags

    def add(self, *tags):
        """Add any number of tags to the Container.
    
        Tags are individual strings that serve to differentiate Containers from
        one another. Sometimes preferable to categories.
    
        :Arguments:
           *tags*
              Tags to add. Must be convertable to strings using the str() builtin.
              May also be a list of tags.

        """
        outtags = list()
        for tag in tags:
            if isinstance(tag, list):
                outtags.extend(tag)
            else:
                outtags.append(tag)
        self._containerfile.add_tags(*outtags)
    
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
        return "<Categories({})>".format(self.dict())

    def __str__(self):
        categories = self.dict()
        agg = "Categories"
        majsep = "="
        seplength = len(agg)

        if not categories:
            out = "No Categories"
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for key in categories.keys():
                out = out + "'{}': '{}'\n".format(key, categories[key])
        return out

    def __getitem__(self, key):
        """Get value at given key.

        :Arguments:
            *key*
                key of value to return

        :Returns:
            *value*
                value corresponding to given key
        """
        categories = self._containerfile.get_categories()
        return categories[key]

    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        outdict = {key: value}
        self._containerfile.add_categories(**outdict)

    def __delitem__(self, category):
        """Remove category from Container.
    
        """
        self._containerfile.del_categories(category)

    def __iter__(self):
        return self._containerfile.get_categories().__iter__()

    def __len__(self):
        return len(self._containerfile.get_categories())

    def dict(self):
        """Get all categories for the Container as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories 
        """
        return self._containerfile.get_categories()

    def add(self, *categorydicts, **categories):
        """Add any number of categories to the Container.

        Categories are key-value pairs of strings that serve to differentiate
        Containers from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will replace
        the value for that category.

        :Keywords:
            *categorydict*
                dict of categories to add; keys used as keys, values used as
                values. Both keys and values must be convertible to strings
                using the str() builtin.
            *categories*
                Categories to add. Keyword used as key, value used as value. Both
                must be convertible to strings using the str() builtin.
        
        """
        outcats = dict()
        for categorydict in categorydicts:
            if isinstance(categorydict, dict):
                outcats.update(categorydict)

        outcats.update(categories)
        self._containerfile.add_categories(**outcats)
    
    def remove(self, *categories, **kwargs):
        """Remove categories from Container.
        
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
    
    def keys(self):
        """Get category keys.
    
        :Returns:
            *keys*
                keys present among categories
        """
        return self._containerfile.get_categories().keys()

    def values(self):
        """Get category values.
    
        :Returns:
            *values*
                values present among categories
        """
        return self._containerfile.get_categories().values()

class Universes(Aggregator):
    """Interface to universes.

    """
    def __repr__(self):
        return "<Universes({})>".format(self.list())

    def __str__(self):
        universes = self.list()
        agg = "Universes"
        majsep = "="
        seplength = len(agg)

        if not universes:
            out = "No universes"
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for universe in universes:
                out = out + "'{}'".format(universe)
                if self._containerfile.get_default() == universe:
                    out = out + ' (default)'
                if self._container._uname == universe:
                    out = out + ' (active)'
                out = out + '\n'
        return out

    def __getitem__(self, handle):
        """Attach universe and return a reference to it.

        :Arguments:
            *handle*
                given name for selecting the universe
        
        :Returns:
            *universe*
                a reference to the newly attached universe
        """
        self.activate(handle)

        return self._container.universe

    def add(self, handle, topology, *trajectory):
        """Add a universe definition to the Sim object.

        A universe is an MDAnalysis object that gives access to the details
        of a simulation trajectory. A Sim object can contain multiple universe
        definitions (topology and trajectory pairs), since it is often
        convenient to have different post-processed versions of the same
        raw trajectory.

        Using an existing universe handle will replace the topology and trajectory
        for that definition; selections for that universe will be retained.

        If there is no current default universe, then the added universe will
        become the default.

        :Arguments:
            *handle*
                given name for selecting the universe
            *topology*
                path to the topology file
            *trajectory*
                path to the trajectory file; multiple files may be given
                and these will be used in order as frames for the trajectory

        """
        outtraj = []
        for traj in trajectory:
            if isinstance(traj, list):
                for t in traj:
                    outtraj.append(t)
            else:
                outtraj.append(traj)

        self._containerfile.add_universe(handle, topology, *outtraj)

        if not self.default():
            self.default(handle)
    
    def remove(self, *handle):
        """Remove a universe definition.
    
        Also removes any selections associated with the universe.

        :Arguments:
            *handle*
                name of universe(s) to delete
        """
        for item in handle:
            self._containerfile.del_universe(item)
    
            if self._container._uname == item:
                self._container._universe = None
                self._container._uname = None

            if self.default() == item:
                self._containerfile.update_default()
    
    def list(self):
        """Get handles for all universe definitions as a list.
    
        :Returns:
            *handles*
                list of all universe handles
        """
        return self._containerfile.list_universes()

    def activate(self, handle=None):
        """Make the selected universe active.
        
        Only one universe definition can be active in a Sim at one time. The
        active universe can be accessed from ``Sim.universe``. Stored
        selections for the active universe can be accessed as items in
        ``Sim.selections``.

        If no handle given, the default universe is loaded.

        If a resnum definition exists for the universe, it is applied.
    
        :Arguments:
            *handle*
                given name for selecting the universe; if ``None``, default
                universe selected
        """
        if not handle:
            handle = self._containerfile.get_default()
    
        if handle:
            udef = self._containerfile.get_universe(handle)
            self._container._universe = Universe(udef[0], *udef[1])
            self._container._uname = handle
            self._apply_resnums()
    
    def current(self):
        """Return the name of the currently active universe.

        :Returns:
            *handle*
                name of currently active universe
        """
        return self._container._uname
                
    def deactivate(self):
        """Deactivate the current universe.
        
        Deactivating the current universe may be necessary to conserve
        memory, since the universe can then be garbage collected.

        """
        self._container._universe = None
        self._container._uname = None
    
    def _apply_resnums(self):
        """Apply resnum definition to active universe.
    
        """
        resnums = self._containerfile.get_resnums(self._container._uname)

        if resnums:
            self._container._universe.residues.set_resnum(resnums)

    def resnums(self, handle, resnums):
        """Define resnums for the given universe.

        Resnums are useful for referring to residues by their canonical resid,
        for instance that stored in the PDB. By giving a resnum definition
        for the universe, this definition will be applied to the universe
        on activation.

        Will overwrite existing resnum definition if it exists.

        :Arguments:
            *handle*
                name of universe to apply resnums to
            *resnums*
                list giving the resnum for each residue in the topology, in
                atom index order; giving ``None`` will delete resnum definition
        """
        if resnums is None:
            self._containerfile.del_resnums(handle)
            
        self._containerfile.update_resnums(handle, resnums)
    
        if handle == self._container._uname:
            self._apply_resnums()
        
    def default(self, handle=None):
        """Mark the selected universe as the default, or get the default universe.

        The default universe is loaded on calls to ``Sim.universe`` or
        ``Sim.selections`` when no other universe is attached.

        If no handle given, returns the current default universe.

        :Arguments:
            *handle*
                given name for selecting the universe; if ``None``, default
                universe is unchanged

        :Returns:
            *default*
                handle of the default universe
        """
        if handle:
            self._containerfile.update_default(handle)
        
        return self._containerfile.get_default()
    
    def define(self, handle):
        """Get the topology and trajectory used for the specified universe.

        :Arguments:
            *handle*
                name of universe to get definition for
        
        :Returns:
            *topology*
                path to the topology file
            *trajectory*
                list of paths to trajectory files
        """
        return self._containerfile.get_universe(handle)

class Selections(Aggregator):
    """Selection manager for Sims.

    Selections are accessible as items using their handles. Each time they are
    called, they are regenerated from the universe that is currently active. In this
    way, changes in the universe topology are reflected in the selections.

    """
    def __repr__(self):
        return "<Selections({})>".format({x: self.define(x) for x in self.keys()})

    def __str__(self):
        selections = self.keys()
        agg = "Selections"
        majsep = "="
        minsep = "-"
        subsep = "| "
        seplength = len(agg)

        if not self._container._uname:
            out = "No universe attached; no Selections to show"
        elif not selections:
            out = "No selections for universe '{}'".format(self._container._uname)
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for selection in selections:
                out = out + "'{}'\n".format(selection)
                for item in self.define(selection):
                    out = out + subsep + "'{}'\n".format(item)
                out = out + minsep*seplength + '\n'

        return out

    def __getitem__(self, handle):
        """Get selection as an AtomGroup for given handle and the active universe.
    
        :Arguments:
            *handle*
                name of selection to return as an AtomGroup

        :Returns:
            *AtomGroup*
                the named selection as an AtomGroup of the active universe
                
        """
        return self.asAtomGroup(handle)

    def __setitem__(self, handle, selection):
        """Selection for the given handle and the active universe.
    
        """
        if isinstance(selection, basestring):
            selection = [selection]
        self._containerfile.add_selection(self._container._uname, handle, *selection)
    
    def __iter__(self):
        return self._containerfile.list_selections(self._container._uname).__iter__()

    def __delitem__(self, handle):
        """Remove stored selection for given handle and the active universe.
    
        """
        self._containerfile.del_selection(self._container._uname, handle)

    def add(self, handle, *selection):
        """Add an atom selection for the attached universe.
        
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
        
    def remove(self, *handle):
        """Remove an atom selection for the attached universe.
        
        :Arguments:
            *handle*
                name of selection(s) to remove
        """
        for item in handle:
            self._containerfile.del_selection(self._container._uname, item)
    
    def keys(self):
        """Return a list of all selection handles.
    
        """
        if self._container._uname:
            return self._containerfile.list_selections(self._container._uname)

    def asAtomGroup(self, handle):
        """Get AtomGroup from active universe corresponding to the given named selection.

        :Arguments:
            *handle*
                name of selection to return as an AtomGroup

        :Returns:
            *AtomGroup*
                the named selection as an AtomGroup of the active universe
        
        """
        selstring = self._containerfile.get_selection(self._container._uname, handle)
        return self._container.universe.selectAtoms(*selstring)

    def define(self, handle):
        """Get selection definition for given handle and the active universe.
    
        :Arguments:
            *handle*
                name of selection to get definition of

        :Returns:
            *definition*
                list of strings defining the atom selection
        """
        return self._containerfile.get_selection(self._container._uname, handle)
    
    def copy(self, universe):
        """Copy defined selections of another universe to the active universe.
    
        :Arguments:
            *universe*
                name of universe definition to copy selections from
        """
        if self._container._uname:
            selections = self._containerfile.list_selections(universe)
    
            for sel in selections:
                seldef = self._containerfile.get_selection(universe, sel)
                self._containerfile.add_selection(self._container._uname, sel, *seldef)

class Members(Aggregator):
    """Member manager for Groups.

    """
    def __repr__(self):
        return "<Members({})>".format(self.names())

    def __str__(self):
        names = self.names()
        containertypes = self.containertypes()
        agg = "Members"
        majsep = "="
        seplength = len(agg)

        if not names:
            out = "No Members"
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for i, name, containertype in zip(xrange(len(names)), names, containertypes):
                out = out + "{}\t{}:\t{}\n".format(i, containertype, name)

        return out

    def __getitem__(self, index):
        """Get member corresponding to the given index.
        
        """
        uuids = self._containerfile.get_members_uuid()
        uuid = uuids[index] 

        if isinstance(uuid, basestring):
            try:
                member = self._container._cache[uuid]
            except KeyError:
                memberdet = self._containerfile.get_member(uuid)
                memberpath = os.path.join(memberdet['abspath'], memberdet['name'])
                if memberdet['containertype'] == 'Sim':
                    member = mds.Sim(memberpath)
                elif memberdet['containertype'] == 'Group':
                    member = mds.Group(memberpath)
                self._container._cache[uuid] = member
        elif isinstance(uuid, list):
            member = list()
            for item in uuid:
                try:
                    member.append(self._container._cache[item])
                except KeyError:
                    memberdet = self._containerfile.get_member(item)
                    memberpath = os.path.join(memberdet['abspath'], memberdet['name'])
                    if memberdet['containertype'] == 'Sim':
                        new = mds.Sim(memberpath)
                    elif memberdet['containertype'] == 'Group':
                        new = mds.Group(memberpath)
                    member.append(new)
                    self._container._cache[item] = new

        return member

    def list(self):
        """Return a list of members.

        Note: modifications of this list won't modify the members of the Group!

        """
        #TODO: need to route these through Finder.
        #TODO: need uuid checks
        members = list()
        for uuid in self._containerfile.get_members_uuid():
            try:
                members.append(self._container._cache[uuid])
            except KeyError:
                member = self._containerfile.get_member(uuid)
                memberpath = os.path.join(member['abspath'], member['name'])
                if member['containertype'] == 'Sim':
                    new = mds.Sim(memberpath)
                elif member['containertype'] == 'Group':
                    new = mds.Group(memberpath)
                members.append(new)
                self._container._cache[uuid] = new
        
        return members

    def containertypes(self):
        """Return a list of member containertypes.

        """
        return self._containerfile.get_members_containertype()

    def names(self):
        """Return a list of member names.

        """
        return self._containerfile.get_members_name()

    def add(self, *containers):
        """Add any number of members to the Group.

        :Arguments:
            *containers*
                Sims and/or Groups to be added; may be a list of Sims and/or
                Groups; Sims or Groups can be given as either objects or paths
                to directories that contain object statefiles
        """
        outconts = list()
        for container in containers:
            if isinstance(container, list):
                self.add(*container)
            elif isinstance(container, mds.Sim) or isinstance(container, mds.Group):
                outconts.append(container)
            elif os.path.exists(container):
                cont = filesystem.path2container(container)
                for c in cont:
                    outconts.append(c)

        for container in outconts:
            self._containerfile.add_member(container.uuid, container.name, container.containertype, container.location)
    
    def remove(self, *indices, **kwargs): 
        """Remove any number of members from the Group.
    
        :Arguments:
            *indices*
                the indices of the members to remove

        :Keywords:
            *all*
                When True, remove all members [``False``]

        """
        uuids = self._containerfile.get_members_uuid()
        if not kwargs.pop('all', False):
            uuids = [ uuids[x] for x in indices ]

        self._containerfile.del_member(*uuids)

class Data(Aggregator):
    """Interface to stored data.

    """
    def __repr__(self):
        return "<Data({})>".format(self.list())

    def _repr_html_(self):
        data = self.list()
        agg = "Data"
        if not data:
            out = "No Data"
        else:
            out = "<h3>{}</h3>".format(agg)
            out = out + "<ul style='list-style-type:none'>"
            for datum in data:
                out = out + "<li>{}</li>".format(datum)
            out = out + "</ul>"
        return out

    def __str__(self):
        data = self.list()
        agg = "Data"
        majsep = "="
        seplength = len(agg)

        if not data:
            out = "No Data"
        else:
            out = agg +'\n'
            out = out + majsep*seplength + '\n'
            for datum in data:
                out = out + "'{}'\n".format(datum)
        return out

    def __iter__(self):
        return self.list().__iter__()

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
                datafile path; None if does not exist
            *datafiletype*
                datafile type; either ``persistence.pddatafile`` or 
                ``persistence.npdatafile``

        """
        datafile = None
        datafiletype = None
        for dfiletype in (persistence.pddatafile, persistence.npdatafile, persistence.pydatafile):
            dfile = os.path.join(self._containerfile.get_location(), 
                                    handle, dfiletype) 
            if os.path.exists(dfile):
                datafile = dfile
                datafiletype = dfiletype

        return (datafile, datafiletype)

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
            filename, filetype = self._get_datafile(handle)

            if filename:
                self._datafile = persistence.DataFile(os.path.join(self._containerfile.get_location(), handle), 
                        logger=self._logger, datafiletype=filetype) 
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

            self._makedirs(dirname)
            self._datafile = persistence.DataFile(dirname, logger=self._logger)

            try:
                out = func(self, handle, *args, **kwargs)
            finally:
                del self._datafile

            return out

        return inner

    def __getitem__(self, handle):
        """Get dataset corresponding to given handle(s).

        If dataset doesn't exist, ``None`` is returned.
        
        :Arguments:
            *handle*
                name of data to retrieve; may also be a list of names

        :Returns:
            *data*
                stored data; if *handle* was a list, will be a list
                of equal length with the stored data as members; will yield
                ``None`` if requested data is nonexistent
                
        """
        if isinstance(handle, list):
            out = list()
            for item in handle:
                out.append(self.retrieve(item))
        elif isinstance(handle, basestring):
            out = self.retrieve(handle)

        return out

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

        A data instance can be a pandas object (Series, DataFrame, Panel),
        a numpy array, or a pickleable python object. If the dataset doesn't
        exist, it is added. If a dataset already exists for the given handle,
        it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data structure to store

        """
        self._datafile.add_data('main', data)
        
    def remove(self, handle, **kwargs):
        """Remove a dataset, or some subset of a dataset.

        Note: in the case the whole dataset is removed, the directory
        containing the dataset file (``Data.h5``) will NOT be removed if it
        still contains file(s) after the removal of the dataset file.

        For pandas objects (Series, DataFrame, or Panel) subsets of the whole
        dataset can be removed using keywords such as *start* and *stop* for
        ranges of rows, and *columns* for selected columns.

        :Arguments:
            *handle*
                name of dataset to delete

        :Keywords:
            *where*
                conditions for what rows/columns to remove
            *start* 
                row number to start selection
            *stop*  
                row number to stop selection
            *columns*
                columns to remove
    
        """
        datafile, datafiletype = self._get_datafile(handle)

        if kwargs and datafiletype == persistence.pddatafile:
            self._delete_data(handle, **kwargs)
        elif datafile:
            os.remove(datafile)
            top = self._containerfile.get_location()
            directory = os.path.dirname(datafile)
            while directory != top:
                try: 
                    os.rmdir(directory)
                    directory = os.path.dirname(directory)
                except OSError:
                    break
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

        :Keywords:
            *where*
                conditions for what rows/columns to remove
            *start* 
                row number to start selection
            *stop*  
                row number to stop selection
            *columns*
                columns to remove

        """
        # only called for pandas objects at the moment
        filename, filetype = self._get_datafile(handle)
        self._datafile.datafiletype = filetype
        try:
            self._datafile.del_data('main', **kwargs)
        except NotImplementedError:
            self._logger.info("Dataset '{}' empty; removing.".format(handle))
            datafile = self._get_datafile(handle)[0]

            os.remove(datafile)
            top = self._containerfile.get_location()
            directory = os.path.dirname(datafile)
            while directory != top:
                try: 
                    os.rmdir(directory)
                    directory = os.path.dirname(directory)
                except OSError:
                    break

    @_read_datafile
    def retrieve(self, handle, **kwargs):
        """Retrieve stored data.

        The stored data structure is read from disk and returned. 

        If dataset doesn't exist, ``None`` is returned.

        For pandas objects (Series, DataFrame, or Panel) subsets of the whole
        dataset can be returned using keywords such as *start* and *stop* for
        ranges of rows, and *columns* for selected columns.

        Also for pandas objects, the *where* keyword takes a string as input
        and can be used to filter out rows and columns without loading the full
        object into memory. For example, given a DataFrame with handle 'mydata'
        with columns (A, B, C, D), one could return all rows for columns A and
        C for which column D is greater than .3 with::

            retrieve('mydata', where='columns=[A,C] & D > .3')

        Or, if we wanted all rows with index = 3 (there could be more than
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
    
        The object must be of the same pandas class (Series, DataFrame, Panel)
        as the existing dataset, and it must have exactly the same columns
        (names included).

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
        top = self._containerfile.get_location()
        for root, dirs, files in os.walk(top):
            if (persistence.pddatafile in files) or (persistence.npdatafile in files) or (persistence.pydatafile in files):
                datasets.append(os.path.relpath(root, start=top))

        datasets.sort()
    
        return datasets

    def locate(self, handle):
        """Get directory location for a stored dataset.

        :Arguments:
            *handle*
                name of data to retrieve location of

        :Returns:
            *datadir*
                absolute path to directory containing stored data

        """
        return os.path.dirname(self._get_datafile(handle)[0])

    def make_filepath(self, handle, filename):
        """Return a full path for a file stored in a data directory, whether
        the file exists or not.

        This is useful if preparing plots or other files derived from the
        dataset, since these can be stored with the data in its own directory.
        This method does the small but annoying work of generating a full path
        for the file.

        This method doesn't care whether or not the path exists; it simply returns
        the path it's asked to build.

        :Arguments:
            *handle*
                name of dataset file corresponds to
            *filename*
                filename of file

        :Returns:
            *filepath*
                absolute path for file

        """
        return os.path.join(os.path.dirname(self._get_datafile(handle)[0]), filename)

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

