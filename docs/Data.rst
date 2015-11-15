==========================
Storing arbitrary datasets
==========================
Treant state files are mainly built to store metadata, but what about storing
(potentially large and time consuming to produce) datasets? Using our Treant
``sprout`` as the example here, say we have generated a `numpy
<http://www.numpy.org/>`__ array of dimension (10^6, 3) that we wish to have
easy access to later ::

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
For other data structures, the Treant will pickle them if it can. The
:class:`datreant.limbs.Data` interface used here is built to make
storage of ~90% of data structures one works with as easy as possible.

Datasets can be nested however you like. For example, say we had several
pandas DataFrames each giving a table of observations for a different subject
from the study the Treant corresponds to. We could just as well make it clear
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
The :class:`datreant.limbs.Data` interface takes advantage of the fact that
Treants are directory trees, giving individual datasets their own place in the
filesystem instead of shoving them into a single file on disk. This is by
design, as it generally means better performance since this means less waiting
for file locks to release from other instances of the same Treant from other
Python sessions. Also, it gives a space to put other files related to the
dataset itself, such as figures produced from it. 

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

Reference: Data
===============
The class :class:`datreant.limbs.Data` is the interface used
by Treants to access their stored datasets. It is not intended to be used
on its own, but is shown here to give a detailed view of its methods.

.. autoclass:: datreant.limbs.Data
    :members:
    :inherited-members:
