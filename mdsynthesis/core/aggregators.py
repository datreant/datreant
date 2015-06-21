"""
Aggregators are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Sim) in
ways that are user-friendly.

In short, an Aggregator is designed to be user friendly on its own, but it can
be used as a backend by a Container, too.

"""
from MDAnalysis import Universe
import os
import numpy as np
from functools import wraps

import persistence
import filesystem
import bundle
import pandas as pd
import mdsynthesis as mds


class Aggregator(object):
    """Core functionality for information aggregators.

    """

    def __init__(self, container, containerfile, logger):
        self._container = container
        self._backend = containerfile
        self._logger = logger

        self._placeholders()

    def _placeholders(self):
        """Initialize any hidden elements.

        """
        pass


class Tags(Aggregator):
    """Interface to tags.

    """

    def __repr__(self):
        return "<Tags({})>".format(self._list())

    def __str__(self):
        tags = self._list()
        agg = "Tags"
        majsep = "="
        seplength = len(agg)

        if not tags:
            out = "No Tags"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for i in xrange(len(tags)):
                out = out + "'{}'\n".format(tags[i])
        return out

    def __iter__(self):
        return self._backend.get_tags().__iter__()

    def __len__(self):
        return len(self._backend.get_tags())

    def _list(self):
        """Get all tags for the Container as a list.

        :Returns:
            *tags*
                list of all tags
        """
        tags = self._backend.get_tags()
        tags.sort()
        return tags

    def add(self, *tags):
        """Add any number of tags to the Container.

        Tags are individual strings that serve to differentiate Containers from
        one another. Sometimes preferable to categories.

        :Arguments:
           *tags*
              Tags to add. Must be convertable to strings using the str()
              builtin.  May also be a list of tags.

        """
        outtags = list()
        for tag in tags:
            if isinstance(tag, list):
                outtags.extend(tag)
            else:
                outtags.append(tag)
        self._backend.add_tags(*outtags)

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
        self._backend.del_tags(*tags, **kwargs)


class Categories(Aggregator):
    """Interface to categories.

    """

    def __repr__(self):
        return "<Categories({})>".format(self._dict())

    def __str__(self):
        categories = self._dict()
        agg = "Categories"
        majsep = "="
        seplength = len(agg)

        if not categories:
            out = "No Categories"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
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
        categories = self._backend.get_categories()
        return categories[key]

    def __setitem__(self, key, value):
        """Set value at given key.

        :Arguments:
            *key*
                key of value to set
        """
        outdict = {key: value}
        self._backend.add_categories(**outdict)

    def __delitem__(self, category):
        """Remove category from Container.

        """
        self._backend.del_categories(category)

    def __iter__(self):
        return self._backend.get_categories().__iter__()

    def __len__(self):
        return len(self._backend.get_categories())

    def _dict(self):
        """Get all categories for the Container as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories

        """
        return self._backend.get_categories()

    def add(self, *categorydicts, **categories):
        """Add any number of categories to the Container.

        Categories are key-value pairs of strings that serve to differentiate
        Containers from one another. Sometimes preferable to tags.

        If a given category already exists (same key), the value given will
        replace the value for that category.

        :Keywords:
            *categorydict*
                dict of categories to add; keys used as keys, values used as
                values. Both keys and values must be convertible to strings
                using the str() builtin.
            *categories*
                Categories to add. Keyword used as key, value used as value.
                Both must be convertible to strings using the str() builtin.

        """
        outcats = dict()
        for categorydict in categorydicts:
            if isinstance(categorydict, dict):
                outcats.update(categorydict)
            else:
                raise TypeError("Invalid arguments; non-keyword" +
                                " arguments must be dicts")

        outcats.update(categories)
        self._backend.add_categories(**outcats)

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
        self._backend.del_categories(*categories, **kwargs)

    def keys(self):
        """Get category keys.

        :Returns:
            *keys*
                keys present among categories
        """
        return self._backend.get_categories().keys()

    def values(self):
        """Get category values.

        :Returns:
            *values*
                values present among categories
        """
        return self._backend.get_categories().values()


class Universes(Aggregator):
    """Interface to universes.

    """

    def __repr__(self):
        return "<Universes({})>".format(self._list())

    def __str__(self):
        universes = self.list()
        agg = "Universes"
        majsep = "="
        seplength = len(agg)

        if not universes:
            out = "No universes"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for universe in universes:
                out = out + "'{}'".format(universe)
                if self._backend.get_default() == universe:
                    out = out + ' (default)'
                if self._container._uname == universe:
                    out = out + ' (active)'
                out = out + '\n'
        return out

    def __contains__(self, item):
        return (item in self._backend.list_universes())

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

        Using an existing universe handle will replace the topology and
        trajectory for that definition; selections for that universe will be
        retained.

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

        self._backend.add_universe(handle, topology, *outtraj)

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
            self._backend.del_universe(item)

            if self._container._uname == item:
                self._container._universe = None
                self._container._uname = None

            if self.default() == item:
                self._backend.update_default()

    def _list(self):
        """Get handles for all universe definitions as a list.

        :Returns:
            *handles*
                list of all universe handles
        """
        return self._backend.list_universes()

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
            handle = self._backend.get_default()

        if handle:
            uh = filesystem.Universehound(self, handle)
            paths = uh.fetch()
            topology = paths['top'][0]
            trajectory = paths['traj']

            self._container._universe = Universe(topology, *trajectory)
            self._container._uname = handle
            self._apply_resnums()

            # update the universe definition; will automatically build current
            # path variants for each file
            self._backend.add_universe(handle, topology, *trajectory)

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
        resnums = self._backend.get_resnums(self._container._uname)

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
            self._backend.del_resnums(handle)

        self._backend.update_resnums(handle, resnums)

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
            self._backend.update_default(handle)

        return self._backend.get_default()

    def define(self, handle, pathtype='abspath'):
        """Get the stored path to the topology and trajectory used for the
        specified universe.

        *Note*: Does no checking as to whether these paths are valid. To
                check this, try activating the universe.

        :Arguments:
            *handle*
                name of universe to get definition for

        :Keywords:
            *pathtype*
                type of path to return; 'abspath' gives an absolute path,
                'relCont' gives a path relative to the Sim's state file

        :Returns:
            *topology*
                path to the topology file
            *trajectory*
                list of paths to trajectory files
        """
        topology, trajectory = self._backend.get_universe(handle)

        return topology[pathtype][0], trajectory[pathtype].tolist()


class Selections(Aggregator):
    """Selection manager for Sims.

    Selections are accessible as items using their handles. Each time they are
    called, they are regenerated from the universe that is currently active. In
    this way, changes in the universe topology are reflected in the selections.

    """

    def __repr__(self):
        return "<Selections({})>".format(
                {x: self.define(x) for x in self.keys()})

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
            out = "No selections for universe '{}'".format(
                self._container._uname)
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for selection in selections:
                out = out + "'{}'\n".format(selection)
                for item in self.define(selection):
                    out = out + subsep + "'{}'\n".format(item)
                out = out + minsep * seplength + '\n'

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
        self._backend.add_selection(
            self._container._uname, handle, *selection)

    def __iter__(self):
        return self._backend.list_selections(
                self._container._uname).__iter__()

    def __delitem__(self, handle):
        """Remove stored selection for given handle and the active universe.

        """
        try:
            self._backend.del_selection(self._container._uname, handle)
        except KeyError:
            raise KeyError(
                    "No such selection '{}'; add it first.".format(handle))

    def add(self, handle, *selection):
        """Add an atom selection for the attached universe.

        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        If a selection with the given *handle* already exists, it is replaced.

        :Arguments:
            *handle*
                name to use for the selection
            *selection*
                selection string; multiple strings may be given and their
                order will be preserved, which is useful for e.g. structural
                alignments
        """
        self._backend.add_selection(
            self._container._uname, handle, *selection)

    def remove(self, *handle):
        """Remove an atom selection for the attached universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection(s) to remove
        """
        for item in handle:
            try:
                self._backend.del_selection(self._container._uname, item)
            except KeyError:
                raise KeyError(
                        "No such selection '{}'; add it first.".format(item))

    def keys(self):
        """Return a list of all selection handles.

        """
        if self._container._uname:
            return self._backend.list_selections(self._container._uname)

    def asAtomGroup(self, handle):
        """Get AtomGroup from active universe from the given named selection.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection to return as an AtomGroup

        :Returns:
            *AtomGroup*
                the named selection as an AtomGroup of the active universe
        """
        try:
            selstring = self._backend.get_selection(
                self._container._uname, handle)
        except KeyError:
            raise KeyError(
                    "No such selection '{}'; add it first.".format(handle))

        return self._container.universe.selectAtoms(*selstring)

    def define(self, handle):
        """Get selection definition for given handle and the active universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection to get definition of

        :Returns:
            *definition*
                list of strings defining the atom selection
        """
        try:
            selstring = self._backend.get_selection(
                            self._container._uname, handle)
        except KeyError:
            raise KeyError(
                    "No such selection '{}'; add it first.".format(handle))

        return selstring

    def copy(self, universe):
        """Copy defined selections of another universe to the active universe.

        :Arguments:
            *universe*
                name of universe definition to copy selections from
        """
        if self._container._uname:
            try:
                selections = self._backend.list_selections(universe)
            except KeyError:
                raise KeyError("No such universe '{}';".format(universe) +
                               " cannot copy selections.")

            for sel in selections:
                seldef = self._backend.get_selection(universe, sel)
                self._backend.add_selection(
                    self._container._uname, sel, *seldef)


class Members(Aggregator, bundle._CollectionBase):
    """Member manager for Groups.

    """

    def _placeholders(self):
        # member cache
        self._cache = dict()
        self._data = None

    def __repr__(self):
        return "<Members({})>".format(self._list())

    def __str__(self):
        names = self.names
        containertypes = self.containertypes
        agg = "Members"
        majsep = "="
        seplength = len(agg)

        if not names:
            out = "No Members"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for i, name, containertype in zip(xrange(len(names)),
                                              names,
                                              containertypes):
                out = out + "{}\t{}:\t{}\n".format(i, containertype, name)

        return out

    @property
    def data(self):
        """The data of the Container.

        Data are user-generated pandas objects (e.g. Series, DataFrames), numpy
        arrays, or any pickleable python object that are stored in the
        Container for easy recall later.  Each data instance is given its own
        directory in the Container's tree.

        """
        if not self._data:
            self._data = MemberData(self)
        return self._data


class MemberAgg(object):
    """Core functionality for aggregators attached to the Members aggregator.

    """

    def __init__(self, members):
        self._members = members


class MemberData(MemberAgg):
    """Manipulators for member data.

    """
    def _list(self, mode='any'):
        """List available datasets.

        :Arguments:
            *mode*
                'any' returns a list of all handles present in at least one
                member; 'all' returns only handles that are present in all
                members

        :Returns:
            *handles*
                list of handles to available datasets

        """
        datasets = [set(member.data) for member in self._members]
        if mode == 'any':
            out = set.union(*datasets)
        elif mode == 'all':
            out = set.intersection(*datasets)

        return list(out)

    # TODO: needs to work for more than just dataframes, series
    def retrieve(self, handle, **kwargs):
        """Retrieve aggregated dataset from all members.

        The stored data structure for each member is read from disk
        and aggregated. The aggregation scheme is dependent on the
        form of the data structure.

        See :meth:`Data.retrieve` for more information on keyword usage.

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
                aggregated data

        """
        agg = None
        for member in self._members:
            d = member.data.retrieve(handle, **kwargs)
            label = len(d.index)*[member.name]
            index = pd.MultiIndex.from_arrays([label, d.index])
            # FIXME: BROKEN!
            d.index = index

            if agg is not None:
                agg = agg.append(d)
            else:
                agg = d

        return agg


class Data(Aggregator):
    """Interface to stored data.

    """

    def __repr__(self):
        return "<Data({})>".format(self._list())

    def _repr_html_(self):
        data = self._list()
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
        data = self._list()
        agg = "Data"
        majsep = "="
        seplength = len(agg)

        if not data:
            out = "No Data"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for datum in data:
                out = out + "'{}'\n".format(datum)
        return out

    def __iter__(self):
        return self._list().__iter__()

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
        for dfiletype in (persistence.pddatafile, persistence.npdatafile,
                          persistence.pydatafile):
            dfile = os.path.join(self._backend.get_location(),
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
                self._datafile = persistence.DataFile(
                        os.path.join(self._backend.get_location(),
                                     handle),
                        logger=self._logger,
                        datafiletype=filetype)
                try:
                    out = func(self, handle, *args, **kwargs)
                finally:
                    del self._datafile
            else:
                self._logger.warning(
                    "No data named '{}' present.".format(handle))
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
            dirname = os.path.join(self._backend.get_location(), handle)

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
            top = self._backend.get_location()
            directory = os.path.dirname(datafile)
            while directory != top:
                try:
                    os.rmdir(directory)
                    directory = os.path.dirname(directory)
                except OSError:
                    break
        else:
            self._logger.info(
                "No data named '{}' present. Nothing to remove".format(handle))

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
            top = self._backend.get_location()
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

    def _list(self):
        """List available datasets.

        :Returns:
            *handles*
                list of handles to available datasets

        """
        datasets = list()
        top = self._backend.get_location()
        for root, dirs, files in os.walk(top):
            if ((persistence.pddatafile in files) or
                    (persistence.npdatafile in files) or
                    (persistence.pydatafile in files)):
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

        This method doesn't care whether or not the path exists; it simply
        returns the path it's asked to build.

        :Arguments:
            *handle*
                name of dataset file corresponds to
            *filename*
                filename of file

        :Returns:
            *filepath*
                absolute path for file

        """
        return os.path.join(os.path.dirname(self._get_datafile(handle)[0]),
                            filename)
