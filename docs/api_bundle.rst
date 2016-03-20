Treant aggregation
==================
These are the API components of ``datreant.core`` for working with multiple
Treants at once, and treating them in aggregate.

.. _Bundle_api:

Bundle
------
The class :class:`datreant.core.Bundle` functions as an ordered set of
Treants. It allows common operations on Treants to be performed in aggregate,
but also includes mechanisms for filtering and grouping based on Treant
attributes, such as tags and categories.


Bundles can be created from all Treants found in a directory tree with
:func:`datreant.core.discover`:

.. autofunction:: datreant.core.discover

They can also be created directly from any number of Treants:

.. autoclass:: datreant.core.Bundle
    :members:
    :inherited-members:
    :noindex:

AggTags
```````
The class :class:`datreant.core.agglimbs.AggTags` is the interface used by
Bundles to access their members' tags.

.. autoclass:: datreant.core.agglimbs.AggTags
    :members:
    :inherited-members:
    :noindex:

AggCategories
`````````````
The class :class:`datreant.core.agglimbs.AggCategories` is the interface used
by Bundles to access their members' categories.

.. autoclass:: datreant.core.agglimbs.AggCategories
    :members:
    :inherited-members:
    :noindex:
