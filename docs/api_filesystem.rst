Filesystem manipulation
=======================
The components of ``datreant.core`` documented here are those designed for
working directly with filesystem objects, namely directories and files.

.. _Tree_api:

Tree
----
The class :class:`datreant.core.Tree` is an interface to a directory in the
filesystem.

.. autoclass:: datreant.core.Tree
    :members:
    :inherited-members:

.. _Leaf_api:

Leaf
----
The class :class:`datreant.core.Leaf` is an interface to a file in the
filesystem.

.. autoclass:: datreant.core.Leaf
    :members:
    :inherited-members:

.. _View_api:

View
----
The class :class:`datreant.core.View` is an ordered set of Trees and Leaves.
It allows for convenient operations on its members as a collective, as well
as providing mechanisms for filtering and subselection.

.. autoclass:: datreant.core.View
    :members:
    :inherited-members:
