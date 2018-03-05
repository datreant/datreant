Treant aggregation
==================
These are the API components of ``datreant`` for working with multiple
Treants at once, and treating them in aggregate.

.. _Bundle_api:

Bundle
------
The class :class:`datreant.Bundle` functions as an ordered set of
Treants. It allows common operations on Treants to be performed in aggregate,
but also includes mechanisms for filtering and grouping based on Treant
attributes, such as tags and categories.


Bundles can be created from all Treants found in a directory tree with
:func:`datreant.discover`:

.. autofunction:: datreant.discover

They can also be created directly from any number of Treants:

.. autoclass:: datreant.Bundle
    :members:
    :inherited-members:

AggTags
```````
The class :class:`datreant.metadata.AggTags` is the interface used by
Bundles to access their members' tags.

.. autoclass:: datreant.metadata.AggTags
    :members:
    :inherited-members:

AggCategories
`````````````
The class :class:`datreant.metadata.AggCategories` is the interface used
by Bundles to access their members' categories.

.. autoclass:: datreant.metadata.AggCategories
    :members:
    :inherited-members:
