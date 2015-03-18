=======================
Datasets and Containers
=======================
MDSynthesis is not an analysis code. On its own, it does not produce output
data given raw simulation data as input. Its scope is limited to the boring
but tedious task of data management and storage. It is intended to bring
value to analysis results by making them easily accessible now and later.

As such, the basic functionality of MDSynthesis is condensed into only two
objects, sometimes referred to as *Containers* in the documentation. These are
the :doc:`Sim <Sim>` and :doc:`Group <Group>` objects.

In brief, a **Sim** is designed to manage and give access to the data corresponding
to a single simulation (the raw trajectory(s), as well as analysis results); a
**Group** gives access to any number of **Sim** or **Group** objects it has as
members (including perhaps itself), and can store analysis results that pertain
to these members collectively. Both types of Container store their underlying
data persistently to disk on the fly. The file locking needed for each
transaction is handled automatically, so more than one python process can be
working with any number of instances of the same Container at the same time.

.. warning:: File locking is generally process safe, but not thread safe. Don't
             use multithreading and try to modify Container elements at the
             same time. Multiprocessing, however, should work just fine.

Persistence as a feature
========================
Containers store their data as directory structures in the file system. Generating
a new **Sim**, for example, with the following ::
    
    >>> # python session 1
    >>> import mdsynthesis as mds
    >>> s = mds.Sim('marklar')

creates a directory called ``marklar`` in the current working directory. It contains
a single file at the moment ::

    > # shell 
    > ls marklar
    Sim.2b4b5800-48a7-4814-ba6d-1e631a09a199.h5

This is the state file containing all the information needed to regenerate an
identical instance of this **Sim**. In fact, we can open a separate python
session (go ahead!) and regenerate this **Sim** immediately there ::

    >>> # python session 2
    >>> import mdsynthesis as mds
    >>> s = mds.Sim('marklar')

Making a modification to the **Sim** in one session, perhaps by adding a tag,
will be reflected in the **Sim** in the other session ::

    >>> # python session 1
    >>> s.tags.add('TIP4P')

    >>> # python session 2
    >>> s.tags
    <Tags(['TIP4P'])>

This is because both objects pull their identifying information from the same
file on disk; they store almost nothing in memory.

.. note:: The uuid of the **Sim** in this example will certainly differ from
          any **Sims** you generate. This is used to differentiate **Sims**
          from each other. Unexpected and broken behavior will result from
          changing the names of state files!

Storing arbitrary datasets
==========================
More on things like tags later, but we really care about storing (potentially
large and time consuming to produce) datasets. Using our **Sim** ``marklar``
as the example here, say we have generated a numpy array of dimension 
(10^6, 3) that gives the minimum distance between the sidechains of three
residues with those of a fourth for each frame in a trajectory ::

    >>> a.shape
    (1000000, 3)

We can store this easily ::

    >>> s.data.add('distances', a)
    >>> s.data
    <Data(['distances'])>

and recall it ::

    >>> s.data['distances'].shape
    (1000000, 3)

Looking at the contents of the directory ``marklar``, we see it has a new
subdirectory corresponding to the name of our stored dataset ::

    > # shell
    > ls marklar
    distances  Sim.h5

which has its own contents ::

    > ls marklar/distances
    npData.h5

This is the data we stored, serialized to disk in the efficient `HDF5
<http://www.hdfgroup.org/HDF5/>`__ data format. Containers will also
store `pandas <http://pandas.pydata.org/>`__ objects using this format.
For other data structures, the Container will pickle them if it can.

Datasets can be nested however you like. For example, say we had several
pandas **DataFrames** each giving the distance with time of each cation in the
simulation with respect to some residue of interest on our protein. We
could just as well make it clear to ourselves that these are all similar
datasets by grouping them together ::

    >>> s.data.add('cations/residue1', df1)
    >>> s.data.add('cations/residue2', df2)
    >>> # we can also use setitem syntax
    >>> s.data['cations/residue3'] = df3
    >>> s.data
    <Data(['cations/residue1', 'cations/residue2', cations/residue3', 
           'distances'])>

and their locations in the filesystem reflect this structure.

Minimal blobs
=============
Individual datasets get their own place in the filesystem instead of all being
shoved into a single file on disk. This is by design, as it generally means
better performance since this means less waiting for file locks to release from
other Container instances. Also, it gives a space to put other files related to
the dataset itself, such as figures produced from it.

You can get the location on disk of a dataset with ::

    >>> s.data.locate('cations/residue1')
    '/home/bob/marklar/cations/residue1'

which is particularly useful for outputting figures.

Another advantage of organizing Containers at the filesystem level is that
datasets can be handled at the filesystem level. Removing a dataset with a ::

    > # shell
    > rm -r marklar/cations/residue2

is immediately reflected by the Container ::

    >>> s.data
    <Data(['cations/residue1', 'cations/residue3', 'distances'])>
    
Datasets can likewise be moved within the Container's directory tree and they
will still be found, with names matching their location relative to the state
file.

Reference: Data
===============
The class :class:`mdsynthesis.core.aggregators.Data` is the interface used
by Containers to access their stored datasets. It is not intended to be used
on its own, but is shown here to give a detailed view of its methods.

.. autoclass:: mdsynthesis.core.aggregators.Data
    :members:
    :inherited-members:
