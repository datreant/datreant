.. _Bundles:

=================================
Coordinating Treants with Bundles
=================================
Similar to a View, a **Bundle** is an ordered set of Treants that makes it easy
to work with them as a single logical unit. Bundles can be constructed in a
variety of ways, but often with a collection of Treants. If our working directory
has a few Treants in it::

    > ls
    acorns/   beans/   sprout/

We can make a Bundle with ::
    
    >>> import datreant.core as dtr
    >>> b = dtr.Bundle('sprout', 'beans', 'acorns')
    >>> b
    <Bundle([<Treant: 'sprout'>, <Treant: 'beans'>, <Treant: 'acorns'>])>

Bundles can also be initialized from existing Treant instances, in addition to
their paths in the filesystem, so ::

    >>> t = dtr.Treant('sprout')
    >>> b = dtr.Bundle(t, 'beans', 'acorns')

would work equally well.


Gathering Treants from the filesystem
=====================================
It can be tedious manually hunting for existing Treants throughout the
filesystem. For this reason the :func:`datreant.core.discover` function
can do this work for us::

    >>> b = dtr.discover('.')
    >>> b
    <Bundle([<Treant: 'beans'>, <Treant: 'sprout'>, <Treant: 'acorns'>])>

For this simple example all our Treants were in this directory, so it's not
quite as useful. But for a directory structure that is deep and convoluted perhaps
from a project spanning years, ``discover`` lets you get a Bundle of all Treants
in the tree with little effort. You can then filter on tags and categories to
get Bundles of the Treants you actually want to work with (see below).

.. autofunction:: datreant.core.discover


Basic member selection
======================



Filtering on Treant tags
========================
Treants are more useful than plain Trees because they carry distinguishing
characteristics beyond just their path in the filesystem. Tags are one of these
distinguishing features, and Bundles can use them directly to filter their
members.

    



Grouping with Treant categories
===============================


Reference: Bundle
=================
.. autoclass:: datreant.core.Bundle
    :members:
    :inherited-members:
