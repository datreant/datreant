=============================================
Filesystem manipulation with Trees and Leaves
=============================================
A Treant functions as a specially marked directory, having a state file with
identifying information that makes it a Treant at all. What's a Treant without
a state file? It's just a **Tree**.

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

Working with paths
==================

getitem, existence, making


A Treant is a Tree
==================

globbing, drawing, 


Reference: Tree
===============
.. autoclass:: datreant.core.Tree
    :members:
    :inherited-members:

Reference: Leaf
===============
.. autoclass:: datreant.core.Leaf
    :members:
    :inherited-members:
