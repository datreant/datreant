.. _Views:

======================================================
Using Views to work with Trees and Leaves collectively
======================================================
A **View** makes it possible to work with arbitrary Trees and Leaves as a
single logical unit. It is an ordered set of its members.


Building a View and selecting members
=====================================
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

But we can also use fancy indexing (which can be useful for re-ordering
members) ::

    >>> v[[2, 0]]
    <View([<Leaf: 'curly.txt'>, <Tree: 'moe/'>])>

Or boolean indexing ::

    >>> v[[True, False, True]]
    <View([<Tree: 'moe/'>, <Leaf: 'curly.txt'>])>

As well as indexing by name ::

    >>> v['curly.txt']
    <View([<Leaf: 'curly.txt'>])>

Note that since the name of a file or directory need not be unique, this always
returns a View.


A View is an abstract Tree...kind of
====================================
A View 



Views from a Tree
=================

Reference: View
===============
.. autoclass:: datreant.core.View
    :members:
    :inherited-members:
