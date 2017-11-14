==========================
Creating and using Treants
==========================
**datreant** is not an analysis library. Its scope is limited to the boring but
tedious task of data management and storage. It is intended to bring value to
analysis results by making them easily accessible now and later.

The basic functionality of datreant is condensed into one object: the
**Treant**. Named after the `talking trees of D&D lore 
<http://wikipedia.org/wiki/Treant>`__, Treants are persistent objects
that live as directory trees in the filesystem and store their state information
to disk on the fly. The file locking needed for each transaction is handled
automatically, so more than one python process can be working with any number
of instances of the same Treant at the same time.

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
    >>> import datreant.core as dtr
    >>> s = dtr.Treant('sprout')

creates a directory called ``sprout`` in the current working directory. It contains
a single directory at the moment ::

    > # shell 
    > ls -a sprout
    .  ..  .datreant

This ``.datreant`` directory is what makes ``sprout`` a Treant. On its own it
serves as a marker, but as we'll see later it can also contain metadata
elements distinguishing this Treant from others.  

Treants are persistent. In fact, we can open a separate python
session (go ahead!) and use this Treant immediately there ::

    >>> # python session 2
    >>> import datreant.core as dtr
    >>> s = dtr.Treant('sprout')

Making a modification to the Treant in one session, perhaps by adding a tag,
will be reflected in the Treant in the other session ::

    >>> # python session 1
    >>> s.tags.add('elm')

    >>> # python session 2
    >>> s.tags
    <Tags(['elm'])>

This is because both objects pull their identifying information from the same
place on disk; they store almost nothing in memory.


API Reference: Treant
=====================
See the :ref:`Treant_api` API reference for more details.
