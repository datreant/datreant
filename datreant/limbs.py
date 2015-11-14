"""
Limbs are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Treant) in
ways that are user-friendly.

In short, an Limb is designed to be user friendly on its own, but it can
be used as a backend by a Treant, too.

"""
import os
from functools import wraps

import numpy as np
import pandas as pd

from datreant import data
from datreant import filesystem
from datreant import collections


class Limb(object):
    """Core functionality for information limbs.

    """

    def __init__(self, treant, treantfile, logger):
        self._treant = treant
        self._backend = treantfile
        self._logger = logger

        self._placeholders()

    def _placeholders(self):
        """Initialize any hidden elements.

        """
        pass


class Tags(Limb):
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
        """Get all tags for the Treant as a list.

        :Returns:
            *tags*
                list of all tags
        """
        tags = self._backend.get_tags()
        tags.sort()
        return tags

    def add(self, *tags):
        """Add any number of tags to the Treant.

        Tags are individual strings that serve to differentiate Treants from
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
        """Remove tags from Treant.

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


class Categories(Limb):
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
        """Remove category from Treant.

        """
        self._backend.del_categories(category)

    def __iter__(self):
        return self._backend.get_categories().__iter__()

    def __len__(self):
        return len(self._backend.get_categories())

    def _dict(self):
        """Get all categories for the Treant as a dictionary.

        :Returns:
            *categories*
                dictionary of all categories

        """
        return self._backend.get_categories()

    def add(self, *categorydicts, **categories):
        """Add any number of categories to the Treant.

        Categories are key-value pairs of strings that serve to differentiate
        Treants from one another. Sometimes preferable to tags.

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
        """Remove categories from Treant.

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


class Members(Limb, collections._CollectionBase):
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
        treanttypes = self.treanttypes
        agg = "Members"
        majsep = "="
        seplength = len(agg)

        if not names:
            out = "No Members"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for i, name, treanttype in zip(xrange(len(names)),
                                           names,
                                           treanttypes):
                out = out + "{}\t{}:\t{}\n".format(i, treanttype, name)

        return out


class MemberAgg(object):
    """Core functionality for limbs attached to the Members limb.

    """

    def __init__(self, members):
        self._members = members


class MemberData(MemberAgg):
    """Manipulators for member data.

    """
    def __repr__(self):
        return "<Data({})>".format(self.keys(mode='any'))

    def _repr_html_(self):
        data = self.keys(mode='any')
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

    def __getitem__(self, handle):
        """Retrieve aggreggated dataset from all members.

        Returns datasets indexed according to member uuids.
        See :meth:`MemberData.retrieve` for more information.

        Raises :exc:`KeyError` if dataset doesn't exist for any members.

        :Arguments:
            *handle*
                name of data to retrieve; may also be a list of names

        :Returns:
            *data*
                aggregated data, indexed by member name; if *handle* was a
                list, will be a list of equal length with the aggregated
                datasets as members

        """
        if isinstance(handle, list):
            out = list()
            for item in handle:
                out.append(self.retrieve(item, by='uuid'))
        elif isinstance(handle, basestring):
            out = self.retrieve(handle, by='uuid')

        return out

    def keys(self, mode='any'):
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

        out = list(out)
        out.sort()

        return out

    def retrieve(self, handle, by='uuid', **kwargs):
        """Retrieve aggregated dataset from all members.

        This is a convenience method. The stored data structure for each member
        is read from disk and aggregated. The aggregation scheme is dependent
        on the form of the data structures pulled from each member:

        pandas DataFrames or Series
            the structures are appended together, with a new level added
            to the index giving the member (see *by*) each set of rows
            came from

        pandas Panel or Panel4D, numpy arrays, pickled python objects
            the structures are returned as a dictionary, with keys giving
            the member (see *by*) and each value giving the corresponding
            data structure

        This method tries to do smart things with the data it reads from each
        member. In particular:
            - members for which there is no data with the given handle are
              skipped
            - the lowest-common-denominator data structure is output; this
              means that if all data structures read are pandas DataFrames,
              then a multi-index DataFrame is returned; if some structures are
              pandas DataFrames, while some are anything else, a dictionary is
              returned

        :Arguments:
            *handle*
                name of data to retrieve

        :Keywords:
            *by*
                top-level index of output data structure; 'name' uses member
                names, 'uuid' uses member uuids; if names are not unique,
                it is better to go with 'uuid' ['uuid']

        See :meth:`Data.retrieve` for more information on keyword usage.

        :Keywords for pandas data structures:
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
                aggregated data structure

        """
        def dict2multiindex(agg):
            agg_mi = None
            for member in agg:
                d = agg[member]
                label = len(d.index)*[member]
                index = pd.MultiIndex.from_arrays([label, d.index])
                d.index = index

                if agg_mi is not None:
                    agg_mi = agg_mi.append(d)
                else:
                    agg_mi = d

            return agg_mi

        # first, check for existence in any member
        if handle not in self.keys('any'):
            raise KeyError(
                    "No dataset '{}' found in any member".format(handle))

        # get indexer from *by* keyword
        if by == 'uuid':
            def get_index(member): return member.uuid
        elif by == 'name':
            def get_index(member): return member.name
            names = [member.name for member in self._members]
            if len(set(names)) != len(names):
                self._members._logger.warning(
                        "Member names not unique; data structure may not" +
                        " look as expected. Set *by* to 'uuid' to avoid this.")
        else:
            raise ValueError(
                    "*by* keyword must be either 'name' or 'uuid'")

        # first, collect all the data into a dictionary, the
        # lowest-common-denominator aggregation structure
        agg = dict()
        for member in self._members:
                agg[get_index(member)] = member.data.retrieve(handle, **kwargs)

        # if data are all Series or all DataFrames, we build a multi-index
        # Series or DataFrame (respectively)
        all_s = all([isinstance(d, pd.Series) for d in agg.values()])
        all_df = all([isinstance(d, pd.DataFrame) for d in agg.values()])

        if all_s or all_df:
            agg = dict2multiindex(agg)

        return agg


class Data(Limb):
    """Interface to stored data.

    """

    def __repr__(self):
        return "<Data({})>".format(self.keys())

    def _repr_html_(self):
        data = self.keys()
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
        data = self.keys()
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
        return self.keys().__iter__()

    def _makedirs(self, p):
        """Make directories and all parents necessary.

        :Arguments:
            *p*
                directory path to make
        """
        try:
            os.makedirs(p)
        except OSError:
            pass

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
        for dfiletype in (data.pddata.pddatafile, data.npdata.npdatafile,
                          data.pydata.pydatafile):
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
        in which data appears, goes missing, etc. while a Treant is loaded.

        ``Note``: methods wrapped with this decorator need to have *handle*
        as the first argument.

        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            filename, filetype = self._get_datafile(handle)

            if filename:
                self._datafile = data.DataFile(
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
        in which data appears, goes missing, etc. while a Treant is loaded.

        ``Note``: methods wrapped with this decorator need to have *handle*
        as the first argument.

        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            dirname = os.path.join(self._backend.get_location(), handle)

            self._makedirs(dirname)
            self._datafile = data.DataFile(dirname, logger=self._logger)

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
        """Store data in Treant.

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

        if kwargs and datafiletype == data.pddata.pddatafile:
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
        side-effects for other instances of this Treant.

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

    def keys(self):
        """List available datasets.

        :Returns:
            *handles*
                list of handles to available datasets

        """
        datasets = list()
        top = self._backend.get_location()
        for root, dirs, files in os.walk(top):
            if ((data.pddata.pddatafile in files) or
                    (data.npdata.npdatafile in files) or
                    (data.pydata.pydatafile in files)):
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
