Filesystem manipulation
=======================
The components of ``datreant`` documented here are those designed for
working directly with filesystem objects, namely directories and files.

.. _Tree_api:

Tree
----
The class :class:`datreant.Tree` is an interface to a directory in the
filesystem.

.. autoclass:: datreant.Tree
    :members:
    :inherited-members:

.. _Leaf_api:

Leaf
----
The class :class:`datreant.Leaf` is an interface to a file in the
filesystem.

.. autoclass:: datreant.Leaf
    :members:
    :inherited-members:

.. _View_api:

View
----
The class :class:`datreant.View` is an ordered set of Trees and Leaves.
It allows for convenient operations on its members as a collective, as well
as providing mechanisms for filtering and subselection.

.. autoclass:: datreant.View
    :members:
    :inherited-members:
