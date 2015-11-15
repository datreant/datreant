==========================
Creating and using Treants
==========================
**datreant** is not an analysis code. Its scope is limited to the boring but
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

What goes into a state file?
============================
The state file of a Treant contains the core pieces of information that define
it. A few of these things are defined in the filesystem itself, including ::

    /home/bob/research/arborea/sprout/Treant.2b4b5800-48a7-4814-ba6d-1e631a09a199.h5
    |_________________________|______|______|____________________________________|__|
              location          name     ^                 uuid                   ^
    |________________________________|   |                                        |
                basedir               treanttype                        statefiletype
    |________________________________________________________________________________|
                                      filepath

This means that changing the location or name of a Treant can be done at the
filesystem level. Although this means that one can change the treanttype and
uuid as well, this is generally not recommended.

Other components, such as the Treant's tags and categories, are stored internally in
the state file (see :ref:`tags_categories` for more on these).
    

Reference: Treant
=================
.. autoclass:: datreant.Treant
    :members:
    :inherited-members:

