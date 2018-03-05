.. _Views:

======================================================
Using Views to work with Trees and Leaves collectively
======================================================
A **View** makes it possible to work with arbitrary Trees and Leaves as a
single logical unit. It is an immutable, ordered set of its members.


.. _Views_selecting:

Building a View and selecting members
=====================================
Views can be built from a list of paths, existing or not. Taking our working
directory with ::

    > ls
    moe/   larry/   curly.txt

We can build a View immediately ::

    >>> import datreant as dtr
    >>> import glob
    >>> v = dtr.View(glob.glob('*'))
    >>> v
    <View(['moe/', 'larry/', 'curly.txt'])>

And we can get to work using it. Since Views are firstly a collection of
members, individual members can be accessed through indexing and slicing ::

    >>> v[1]
    <Tree: 'larry/'>

    >>> v[1:]
    <View(['larry/', 'curly.txt'])>

But we can also use fancy indexing (which can be useful for getting a
re-ordering of members) ::

    >>> v[[2, 0]]
    <View(['curly.txt', 'moe/'])>

Or boolean indexing ::

    >>> v[[True, False, True]]
    <View(['moe/', 'curly.txt'])>

As well as indexing by name ::

    >>> v['curly.txt']
    <View(['curly.txt'])>

Note that since the name of a file or directory need not be unique, this always
returns a View.


Filtering View members
======================
Often we might obtain a View of a set of files and directories and then use
the View itself to filter down into the set of things we actually want. There
are a number of convenient ways to do this.

Want only the Trees? ::

    >>> v.membertrees
    <View(['moe/', 'larry/'])>

Or only the Leaves? ::

    >>> v.memberleaves
    <View(['curly.txt'])>

We can get more granular and filter members using glob patterns on their names::

    >>> v.globfilter('*r*')
    <View(['larry/', 'curly.txt'])>

And since all these properties and methods return Views, we can stack
operations::

    >>> v.globfilter('*r*').memberleaves
    <View(['curly.txt'])>


.. _Views_from_Tree:

Views from a Tree
=================
A common use of a View is to introspect the children of a Tree. If we have a
look inside one of our directories ::

    > ls moe/
    about_moe.txt  more_moe.pdf  sprout/

We find two files and a directory. We can get at the files with ::

    >>> moe = v['moe'][0]
    >>> moe.leaves()
    <View(['about_moe.txt', 'more_moe.pdf'])>


and the directories with ::

    >>> moe.trees()
    <View(['sprout/'])>

Both these properties leave out hidden children by default, since hidden files
are often hidden to keep them out of the way of most work. But we can get at
these easily, too::

    >>> moe.trees(hidden=True)
    <View(['.hiding_here/', 'sprout/'])>

Want all the children? ::

    >>> moe.children(hidden=True)
    <View(['sprout/', 'about_moe.txt', 'more_moe.pdf', '.hiding_here/'])>


A View is an ordered set
========================
Because a View is a set, adding members that are already present results
in no a new View with nothing additional::

    >>> v = v + Tree('moe')
    >>> v
    <View(['moe/', 'larry/', 'curly.txt'])>

But a View does have a sense of order, so we could, for example, meaningfully
get a View with the order of members reversed::

    >>> v[::-1]
    <View(['curly.txt', 'larry/', 'moe/'])>

Because it is functionally a set, operations between Views work as expected.
Making another View with ::

    >>> v2 = dtr.View('moe', 'nonexistent_file.txt')

we can get the union::

    >>> v | v2
    <View(['moe/', 'larry/', 'curly.txt', 'nonexistent_file.txt'])>

the intersection::

    >>> v & v2
    <View(['moe/'])>

differences::

    >>> v - v2
    <View(['larry/', 'curly.txt'])>

    >>> v2 - v
    <View(['nonexistent_file.txt'])>

or the symmetric difference::

    >>> v ^ v2
    <View(['curly.txt', 'larry/', 'nonexistent_file.txt'])>


Collective properties and methods of a View
===========================================
A View is a collection of Trees and Leaves, but it has methods and properties
that mirror those of Trees and Leaves that allow actions on all of its members
in aggregate. For example, we can directly get all directories and files within
each member Tree::
    
    >>> v.children(hidden=True)
    <View(['sprout/', 'about_moe.txt', 'more_moe.pdf', '.hiding_here',
           'about_larry.txt'])>

Or we could get all children that match a glob pattern::

    >>> v.glob('*moe*')
    <View(['about_moe.txt', 'more_moe.pdf'])>

Note that this is the equivalent of doing something like::

    >>> dtr.View([tree.glob(pattern) for tree in v.membertrees])

In this way, a View functions analogously for Trees and Leaves as a Bundle does
for Treants. See :ref:`Bundles` for more on this theme.

API Reference: View
===================
See the :ref:`View_api` API reference for more details.
