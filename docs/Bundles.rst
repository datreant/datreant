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
All the same selection patterns that work for Views (see :ref:`Views_selecting`).
This includes indexing with integers::

    >>> b = dtr.discover()
    >>> b[1]
    <Treant: 'sprout'>

slicing::

    >>> b[1:]
    <Bundle([<Treant: 'sprout'>, <Treant: 'acorns'>])>

fancy indexing:: 

    >>> b[[1, 2, 0]]
    <Bundle([<Treant: 'sprout'>, <Treant: 'acorns'>, <Treant: 'beans'>])>

boolean indexing::

    >>> b[[False, False, True]]
    <Bundle([<Treant: 'acorns'>])>

and indexing by Treant name::

    >>> b['sprout']
    <Bundle([<Treant: 'sprout'>])>

Note that since Treant names need not be unique, indexing by name always yields
a Bundle.


Filtering on Treant tags
========================
Treants are more useful than plain Trees because they carry distinguishing
characteristics beyond just their path in the filesystem. Tags are one of these
distinguishing features, and Bundles can use them directly to filter their
members.

The aggregated tags for all members in a Bundle are accessible via
:attr:`datreant.core.Bundle.tags`. Just calling this property gives a view of
the tags present in every member Treant::

    >>> b.tags
    <AggTags(['plant'])>

But our Treants probably have more than just this one tag. We can get at the
tags represented by at least one Treant in the Bundle with ::

    >>> b.tags.any
    {'for eating',
     'for squirrels',
     'fruitless',
     'not for humans',
     'plant',
     'plentiful',
     'tiny'}

Since tags function as a set, we get back a set. Likewise we have ::

    >>> b.tags.all
    {'plant'}

which we've already seen.


Using tag expressions to select members
---------------------------------------

We can use getitem syntax to query the members of Bundle. For example, giving a
single tag like ::

    >>> b.tags['for eating']
    [True, False, False]

gives us back a list of booleans. This can be used directly on the Bundle as
a boolean index to get back a subselection of its members::

    >>> b[b.tags['for eating']]
    <Bundle([<Treant: 'beans'>])>

We can also provide multiple tags to match more Treants::

    >>> b[b.tags['for eating', 'not for humans']]
    <Bundle([<Treant: 'beans'>, <Treant: 'acorns'>])>

The above is equivalent to giving a tuple of tags to match, as below::

    >>> b[b.tags[('for eating', 'not for humans')]]
    <Bundle([<Treant: 'beans'>, <Treant: 'acorns'>])>

Using a tuple functions as an "or"-ing of the tags given, in which case
the resulting members are those that have at least one of the tags inside
the tuple.

But if we give a list instead, we get::

    >>> b[b.tags[['for eating', 'not for humans']]]
    <Bundle([])>

...something else, in this case nothing. Giving a list functions as an
"and"-ing of the tags given inside, so the above query will only give members
that have both 'for eating' and 'not for humans' as tags. There were none in
this case. 

Lists and tuples can be nested to build complex and/or selections. In addition,
sets can be used to indicate negation ("not")::

    >>> b[b.tags[{'for eating'}]]
    <Bundle([<Treant: 'sprout'>, <Treant: 'acorns'>])>

Putting multiple tags inside a set functions as a negated "and"-ing of the
contents::

    >>> b[b.tags[{'for eating', 'not for humans'}]]
    <Bundle([<Treant: 'beans'>, <Treant: 'sprout'>, <Treant: 'acorns'>])>


Fuzzy matching for tags
-----------------------


Reference: AggTags
------------------
.. autoclass:: datreant.core.agglimbs.AggTags
    :members:
    :inherited-members:



Grouping with Treant categories
===============================

Reference: AggCategories
------------------------
.. autoclass:: datreant.core.agglimbs.AggCategories
    :members:
    :inherited-members:

Reference: Bundle
=================
.. autoclass:: datreant.core.Bundle
    :members:
    :inherited-members:
