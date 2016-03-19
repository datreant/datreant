.. _Views:

======================================================
Using Views to work with Trees and Leaves collectively
======================================================
A **View** makes it possible to work with arbitrary Trees and Leaves as a
single logical unit. It is an ordered set of its members.


.. _Views_selecting:

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


Filtering View members
======================
Often we might obtain a View of a set of files and directories and then use
the View itself to filter down into the set of things we actually want. There
are a number of convenient ways to do this.

Want only the Trees? ::

    >>> v.membertrees
    <View([<Tree: 'moe/'>, <Tree: 'larry/'>])>

Or only the Leaves? ::

    >>> v.memberleaves
    <View([<Leaf: 'curly.txt'>])>

We can get more granular and filter members using glob patterns on their names::

    >>> v.filter('*r*')
    <View([<Tree: 'larry/'>, <Leaf: 'curly.txt'>])>

And since all these properties and methods return Views, we can stack
operations::

    >>> v.filter('*r*').memberleaves
    <View([<Leaf: 'curly.txt'>])>


.. _Views_from_Tree:

Views from a Tree
=================
A common use of a View is to introspect the children of a Tree. If we have a
look inside one of our directories ::

    > ls moe/
    about_moe.txt  more_moe.pdf  sprout/

We find two files and a directory. We can get at the files with ::

    >>> moe = v['moe'][0]
    >>> moe.leaves
    <View([<Leaf: 'moe/about_moe.txt'>, <Leaf: 'moe/more_moe.pdf'>])>


and the directories with ::

    >>> moe.trees
    <View([<Tree: 'moe/sprout/'>])>

Both these properties leave out hidden children, since hidden files are often
hidden to keep them out of the way of most work. But we can get at these
easily, too::

    >>> moe.hidden
    <View([<Leaf: 'moe/.hiding_here'>])>

Want all the children? ::

    >>> moe.children
    <View([<Tree: 'moe/sprout/'>, <Leaf: 'moe/about_moe.txt'>,
           <Leaf: 'moe/more_moe.pdf'>, <Leaf: 'moe/.hiding_here'>])>


A View is an ordered set
========================
Because a View is a set, adding members that are already present results
in no change to the View::

    >>> v.add('moe')
    >>> v
    <View([<Tree: 'moe/'>, <Tree: 'larry/'>, <Leaf: 'curly.txt'>])>

    # addition of Trees/Leaves to a View also returns a View
    >>> v + dtr.Tree('moe')
    <View([<Tree: 'moe/'>, <Tree: 'larry/'>, <Leaf: 'curly.txt'>])>

But a View does have a sense of order, so we could, for example, meaningfully
get a View with the order of members reversed::

    >>> v[::-1]
    <View([<Leaf: 'curly.txt'>, <Tree: 'larry/'>, <Tree: 'moe/'>])>

Because it is functionally a set, operations between Views work as expected.
Making another View with ::

    >>> v2 = dtr.View('moe', 'nonexistent_file.txt')

we can get the union::

    >>> v | v2
    <View([<Tree: 'moe/'>, <Tree: 'larry/'>, <Leaf: 'curly.txt'>,
           <Leaf: 'nonexistent_file.txt'>])>

the intersection::

    >>> v & v2
    <View([<Tree: 'moe/'>])>

differences::

    >>> v - v2
    <View([<Tree: 'larry/'>, <Leaf: 'curly.txt'>])>

    >>> v2 - v
    <View([<Leaf: 'nonexistent_file.txt'>])>

or the symmetric difference::

    >>> v ^ v2
    <View([<Leaf: 'curly.txt'>, <Tree: 'larry/'>,
           <Leaf: 'nonexistent_file.txt'>])>


Collective properties and methods of a View
===========================================
A View is a collection of Trees and Leaves, but it has methods and properties
that mirror those of Trees and Leaves that allow actions on all of its members
in aggregate. For example, we can directly get all directories and files within
each member Tree::
    
    >>> v.children
    <View([<Tree: 'moe/sprout/'>, <Leaf: 'moe/about_moe.txt'>,
           <Leaf: 'moe/more_moe.pdf'>, <Leaf: 'moe/.hiding_here'>,
           <Leaf: 'larry/about_larry.txt'>])>

Or we could get all children that match a glob pattern::

    >>> v.glob('*moe*')
    <View([<Leaf: 'moe/about_moe.txt'>, <Leaf: 'moe/more_moe.pdf'>])>

Note that this is the equivalent of doing something like::

    >>> dtr.View([tree.glob(pattern) for tree in v.membertrees])

In this way, a View functions analogously for Trees and Leaves as a Bundle does
for Treants. See :ref:`Bundles` for more on this theme.


Reference: View
===============
.. autoclass:: datreant.core.View
    :members:
    :inherited-members:
