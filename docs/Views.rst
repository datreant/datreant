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

.. _Views_from_Tree:

Views from a Tree
=================
A common use of a View is to introspect the children of a Tree. If we have
a look inside one of our directories ::

    > ls moe/
    about_moe.txt  more_moe.pdf  sprout/

We find two files and a directory. We can get at the files with ::

    >>> moe = v['moe'][0]
    >>> moe.childleaves



A View is an abstract Tree...kind of
====================================
A View is roughly duck-typed to behave like an abstract Tree, in which the
contents of all its member Trees are manipulatable as if the View itself
was a Tree with those contents. For example, we can directly get all
directories and files within each Tree that match a glob pattern ::

    >>> 




Reference: View
===============
.. autoclass:: datreant.core.View
    :members:
    :inherited-members:
