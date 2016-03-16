.. _Views:

======================================================
Using Views to work with Trees and Leaves collectively
======================================================
A **View** makes it possible to work with arbitrary Trees and Leaves as a
single logical unit. It is an ordered set of its members.


Building a View
===============
Views can be built from a list of paths, existing or not. Taking our working
directory with ::

    > ls
    moe/   larry/   curly.txt

We can build a View immediately ::

    >>> import datreant.core as dtr
    >>> import glob
    >>> v = dtr.View(glob.glob('*'))
    >>> v
    <View([<Tree: 'moe/'>, <Tree: 'larry/'>, <Leaf: 'curly.txt'>])>

And we can get to work using it. Since Views are firstly a collection of
members, individual members can be accessed through indexing and slicing ::

    >>> v[1]
    <Tree: 'larry/'>

    >>> v[1:]
    <View([<Tree: 'larry/'>, <Leaf: 'curly.txt'>])>



Views from a Tree
=================

Reference: View
===============
.. autoclass:: datreant.core.View
    :members:
    :inherited-members:
