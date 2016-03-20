=============================================
Filesystem manipulation with Trees and Leaves
=============================================
A Treant functions as a specially marked directory, having a state file with
identifying information. What's a Treant without a state file? It's just a
**Tree**.

datreant gives pythonic access to the filesystem by way of **Trees** and
**Leaves** (directories and files, respectively). Say our current working
directory has two directories and a file ::

    > ls
    moe/   larry/   curly.txt

We can use Trees and Leaves directly to manipulate them ::

    >>> import datreant.core as dtr
    >>> t = dtr.Tree('moe')
    >>> t
    <Tree: 'moe'>
    
    >>> l = dtr.Leaf('curly.txt')
    >>> l
    <Leaf: 'curly.txt'>

These objects point to a specific path in the filesystem, which doesn't
necessarily have to exist. Just as with Treants, more than one instance
of a Tree or Leaf can point to the same place.


Working with Trees 
==================
**Tree** objects can be used to introspect downward into their directory
structure. Since a Tree is essentially a container for its own child Trees and
Leaves, we can use getitem syntax to dig around ::

    >>> t = dtr.Tree('moe')
    >>> t['a/directory/']
    <Tree: 'moe/a/directory/'>

    >>> t['a/file']
    <Leaf: 'moe/a/file'>

Paths that resolve as being inside a Tree give `True` for membership tests ::

    >>> t['a/file'] in t
    True

Note that these items need not exist ::

    >>> t['a/file'].exists
    False

in which case whether a Tree or Leaf is returned is dependent on an ending
``/``. We can create directories and empty files easily enough, though::

    >>> adir = t['a/directory/'].make()
    >>> adir.exists
    True

    >>> afile = t['a/file'].make()
    >>> afile.exists
    True

.. note:: For accessing directories and files that exist, getitem syntax isn't
          sensitive to ending ``/`` separators to determine whether to give a
          Tree or a Leaf.

API Reference: Tree
-------------------
See the :ref:`Tree_api` API reference for more details.


A Treant is a Tree
==================
The **Treant** object is a subclass of a Tree, so the above all applies to
Treant behavior. Some methods of Trees are especially useful when working with
Treants. One of these is ``draw`` ::

    >>> s = dtr.Treant('sprout')
    >>> s['a/new/file'].make()
    >>> s['a/.hidden/directory/'].make()
    >>> s.draw()
    sprout/
     +-- Treant.839c7265-5331-4224-a8b6-c365f18b9997.json
     +-- a/
         +-- new/
         |   +-- file
         +-- .hidden/
             +-- directory/

which gives a nice ASCII-fied visual of the Tree. We can also obtain a
collection of Trees and/or Leaves in the Tree with globbing ::

    >>> s.glob('a/*')
    <View([<Tree: 'sprout/a/.hidden/'>, <Tree: 'sprout/a/new/'>])>

See :ref:`Views` for more about the **View** object, and how it can be used to
manipulate many Trees and Leaves as a single logical unit. More details on
how to introspect Trees with Views can be found in :ref:`Views_from_Tree`.


File operations with Leaves
===========================


API Reference: Leaf 
-------------------
See the :ref:`Leaf_api` API reference for more details.

