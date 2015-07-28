====================
Treants and Datasets
====================
datreant is not an analysis code. Its scope is limited to the boring but
tedious task of data management and storage. It is intended to bring value to
analysis results by making them easily accessible now and later.

The basic functionality of datreant is condensed into one object: the
**Treant**. Named after the `talking trees of D&D lore 
<http://wikipedia.org/wiki/Treant>`__, Treants are persistent objects
that live as directory trees in the filesystem. Treants store their
underlying data persistently to disk on the fly. The file locking needed for
each transaction is handled automatically, so more than one python process can
be working with any number of instances of the same Treant at the same
time.

.. warning:: File locking is performed with POSIX advisory locks. These are
             not guaranteed to work perfectly on all platforms and file
             systems, so use caution when changing the stored attributes
             of a Treant in more than one process. Also, though advisory locks
             are mostly process safe, they are definitely not thread safe.
             Don't use multithreading and try to modify Treant elements at the
             same time.

Persistence as a feature
========================
Treants store their data as directory structures in the file system. Generating
a new Treant, for example, with the following ::
    
    >>> # python session 1
    >>> import datreant as dtr
    >>> s = dtr.Treant('sprout')

creates a directory called ``sprout`` in the current working directory. It contains
a single file at the moment ::

    > # shell 
    > ls sprout
    Treant.2b4b5800-48a7-4814-ba6d-1e631a09a199.h5

The name of this file includes the type of Treant it corresponds to, as
well as the ``uuid`` of the Treant, which is its unique identifier. This
is the state file containing all the information needed to regenerate an
identical instance of this Treant. In fact, we can open a separate python
session (go ahead!) and regenerate this Treant immediately there ::

    >>> # python session 2
    >>> import datreant as dtr
    >>> t = dtr.Treant('sprout')

Making a modification to the Treant in one session, perhaps by adding a tag,
will be reflected in the Treant in the other session ::

    >>> # python session 1
    >>> t.tags.add('elm')

    >>> # python session 2
    >>> t.tags
    <Tags(['elm'])>

This is because both objects pull their identifying information from the same
file on disk; they store almost nothing in memory.

.. note:: The ``uuid`` of the Treant in this example will certainly differ from
          any Treants you generate. This is used to differentiate Treants
          from each other. Unexpected and broken behavior will result from
          changing the names of state files!

Storing arbitrary datasets
==========================
More on things like tags later, but we really care about storing (potentially
large and time consuming to produce) datasets. Using our Treant ``sprout``
as the example here, say we have generated a `numpy <http://www.numpy.org/>`__
array of dimension (10^6, 3) that we wish to have easy access to later ::

    >>> a.shape
    (1000000, 3)

We can store this easily ::

    >>> t.data.add('something_wicked', a)
    >>> t.data
    <Data(['something_wicked'])>

and recall it ::

    >>> t.data['something_wicked'].shape
    (1000000, 3)

Looking at the contents of the directory ``sprout``, we see it has a new
subdirectory corresponding to the name of our stored dataset ::

    > # shell
    > ls sprout
    something_wicked  Treant.2b4b5800-48a7-4814-ba6d-1e631a09a199.h5

which has its own contents ::

    > ls sprout/something_wicked
    npData.h5

This is the data we stored, serialized to disk in the efficient `HDF5
<http://www.hdfgroup.org/HDF5/>`__ data format. Treants will also
store `pandas <http://pandas.pydata.org/>`__ objects using this format.
For other data structures, the Treant will pickle them if it can.

Datasets can be nested however you like. For example, say we had several
pandas DataFrames each giving a table of observations for a particular subject
of the study the Treant corresponds to. We could just as well make it clear
to ourselves that these are similar datasets by grouping them together ::

    >>> t.data.add('subjects/leafy', df1)
    >>> t.data.add('subjects/barkley', df2)
    >>> # we can also use setitem syntax
    >>> t.data['cations/twiggy'] = df3
    >>> t.data
    <Data(['subjects/leafy', 'subjects/barkley', subjects/twiggy', 
           'something_wicked'])>

and their locations in the filesystem reflect this structure.

Minimal blobs
=============
Individual datasets get their own place in the filesystem instead of all being
shoved into a single file on disk. This is by design, as it generally means
better performance since this means less waiting for file locks to release from
other Treant instances. Also, it gives a space to put other files related to
the dataset itself, such as figures produced from it.

You can get the location on disk of a dataset with ::

    >>> t.data.locate('subjects/barkley')
    '/home/bob/sprout/subjects/barkley'

which is particularly useful for outputting figures.

Another advantage of organizing Treants at the filesystem level is that
datasets can be handled at the filesystem level. Removing a dataset with a ::

    > # shell
    > rm -r sprout/subjects/leafy

is immediately reflected by the Treant ::

    >>> t.data
    <Data(['subjects/barkley', 'subjects/twiggy', 'something_wicked'])>
    
Datasets can likewise be moved within the Treant's directory tree and they
will still be found, with names matching their location relative to the state
file.

Reference: Treant
==============
.. autoclass:: datreant.Treant
    :members:
    :inherited-members:

Reference: Data
===============
The class :class:`datreant.aggregators.Data` is the interface used
by Treants to access their stored datasets. It is not intended to be used
on its own, but is shown here to give a detailed view of its methods.

.. autoclass:: datreant.aggregators.Data
    :members:
    :inherited-members:
