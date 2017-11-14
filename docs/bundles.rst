.. _Bundles:

=================================
Coordinating Treants with Bundles
=================================
Similar to a View, a **Bundle** is an immutable, ordered set of Treants that
makes it easy to work with them as a single logical unit. Bundles can be
constructed in a variety of ways, but often with a collection of Treants. If
our working directory has a few Treants in it::

    > ls
    elm/  maple/  oak/  sequoia/

We can make a Bundle with ::
    
    >>> import datreant.core as dtr
    >>> b = dtr.Bundle('elm', 'maple', 'oak', 'sequoia')
    >>> b
    <Bundle(['elm', 'maple', 'oak', 'sequoia'])>

Bundles can also be initialized from existing Treant instances, in addition to
their paths in the filesystem, so ::

    >>> t = dtr.Treant('elm')
    >>> b = dtr.Bundle(t, 'maple', 'oak', 'sequoia')

would work equally well.


Gathering Treants from the filesystem
=====================================
It can be tedious manually hunting for existing Treants throughout the
filesystem. For this reason the :func:`~datreant.core.discover` function
can do this work for us::

    >>> b = dtr.discover('.')
    >>> b
    <Bundle(['sequoia', 'maple', 'oak', 'elm'])>

For this simple example all our Treants were in this directory, so it's not
quite as useful. But for a directory structure that is deep and convoluted
perhaps from a project spanning years, :func:`~datreant.core.discover` lets
you get a Bundle of all Treants in the tree with little effort. You can then
filter on tags and categories to get Bundles of the Treants you actually want
to work with.

See the :func:`datreant.core.discover` API reference for more details.

Basic member selection
======================
All the same selection patterns that work for Views (see :ref:`Views_selecting`)
work for Bundles. This includes indexing with integers::

    >>> b = dtr.discover()
    >>> b[1]
    <Treant: 'maple'>

slicing::

    >>> b[1:]
    <Bundle(['maple', 'oak', 'elm'])>

fancy indexing:: 

    >>> b[[1, 2, 0]]
    <Bundle(['maple', 'oak', 'sequoia'])>

boolean indexing::

    >>> b[[False, False, True, False]]
    <Bundle(['oak'])>

and indexing by Treant name::

    >>> b['oak']
    <Bundle(['oak'])>

Note that since Treant names need not be unique, indexing by name always yields
a Bundle.


Filtering on Treant tags
========================
Treants are more useful than plain Trees because they carry distinguishing
characteristics beyond just their path in the filesystem. Tags are one of these
distinguishing features, and Bundles can use them directly to filter their
members.

.. note:: For a refresher on using tags with individual Treants, see 
          :ref:`Tags_guide`. Everything that applies to using tags with
          individual Treants applies to using them in aggregate with Bundles.

The aggregated tags for all members in a Bundle are accessible via
:attr:`datreant.core.Bundle.tags`. Just calling this property gives a view of
the tags present in every member Treant::

    >>> b.tags
    <AggTags(['plant'])>

But our Treants probably have more than just this one tag. We can get at the
tags represented by at least one Treant in the Bundle with ::

    >>> b.tags.any
    {'building',
     'firewood',
     'for building',
     'furniture',
     'huge',
     'paper',
     'plant',
     'shady',
     'syrup'}

Since tags function as a set, we get back a set. Likewise we have ::

    >>> b.tags.all
    {'plant'}

which we've already seen.

Using tag expressions to select members
---------------------------------------
We can use getitem syntax to query the members of Bundle. For example, giving a
single tag like ::

    >>> b.tags['building']
    [False, False, True, True]

gives us back a list of booleans. This can be used directly on the Bundle as
a boolean index to get back a subselection of its members::

    >>> b[b.tags['building']]
    <Bundle(['oak', 'elm'])>

We can also provide multiple tags to match more Treants::

    >>> b[b.tags['building', 'furniture']]
    <Bundle(['maple', 'oak', 'elm'])>

The above is equivalent to giving a tuple of tags to match, as below::

    >>> b[b.tags[('building', 'furniture')]]
    <Bundle(['maple', 'oak', 'elm'])>

Using a tuple functions as an "or"-ing of the tags given, in which case
the resulting members are those that have at least one of the tags inside
the tuple.

But if we give a list instead, we get::

    >>> b[b.tags[['building', 'furniture']]]
    <Bundle([])>

...something else, in this case nothing. Giving a list functions as an
"and"-ing of the tags given inside, so the above query will only give members
that have both 'building' and 'furniture' as tags. There were none in this
case. 

Lists and tuples can be nested to build complex and/or selections. In addition,
sets can be used to indicate negation ("not")::

    >>> b[b.tags[{'furniture'}]]
    <Bundle(['sequoia', 'oak', 'elm'])>

Putting multiple tags inside a set functions as a negated "and"-ing of the
contents::

    >>> b[b.tags[{'building', 'furniture'}]]
    <Bundle(['sequoia', 'maple', 'oak', 'elm'])>

which is the opposite of the empty Bundle we got when we did the "and"-ing of
these tags earlier.

Fuzzy matching for tags
-----------------------
Over the course of a project spanning years, you might add several variations
of essentially the same tag to different Treants. For example, it looks like we
might have two different tags that mean the same thing among the Treants in our
Bundle::

    >>> b.tags
    {'building',
     'firewood',
     'for building',
     'furniture',
     'huge',
     'paper',
     'plant',
     'shady',
     'syrup'}

Chances are good we meant the same thing when we added 'building' and 
'for building' to these Treants. How can we filter on these without explicitly
including each one in a tag expression?

We can use fuzzy matching::

    >>> b.tags.fuzzy('building', scope='any')
    ('building', 'for building')

which we can use directly as an "or"-ing in a tag expression::

    >>> b[b.tags[b.tags.fuzzy('building', scope='any')]]
    <Bundle(['oak', 'elm'])>

The threshold for fuzzy matching can be set with the ``threshold`` parameter.
See the API reference for :meth:`~datreant.core.metadata.AggTags.fuzzy` for more
details on how to use this method.

Grouping with Treant categories
===============================
Besides tags, categories are another mechanism for distinguishing Treants from
each other. We can access these in aggregate with a Bundle, but we can also use
them to build groupings of members by category value.

.. note:: For a refresher on using categories with individual Treants, see 
          :ref:`Categories_guide`. Much of what applies to using categories
          with individual Treants applies to using them in aggregate with
          Bundles.

The aggregated categories for all members in a Bundle are accessible via
:attr:`datreant.core.Bundle.categories`. Just calling this property gives a
view of the categories with keys present in every member Treant::

    >>> b.categories
    <AggCategories({'age': ['adult', 'young', 'young', 'old'], 
                    'type': ['evergreen', 'deciduous', 'deciduous', 'deciduous'], 
                    'bark': ['fibrous', 'smooth', 'mossy', 'mossy']})>

We see that here, the values are lists, with each element of the list giving
the value for each member, in member order. This is how categories behave when
accessing from Bundles, since each member may have a different value for a
given key.

But just as with tags, our Treants probably have more than just the keys 'age',
'type', and 'bark' among their categories. We can get a dictionary of the
categories with each key present among at least one member with ::

    >>> b.categories.any
    {'age': ['adult', 'young', 'young', 'old'],
     'bark': ['fibrous', 'smooth', 'mossy', 'mossy'],
     'health': [None, None, 'good', 'poor'],
     'nickname': ['redwood', None, None, None],
     'type': ['evergreen', 'deciduous', 'deciduous', 'deciduous']}

Note that for members that lack a given key, the value returned in the
corresponding list is ``None``. Since ``None`` is not a valid value for a
category, this unambiguously marks the key as being absent for these members.

Likewise we have ::

    >>> b.categories.all
    {'age': ['adult', 'young', 'young', 'old'],
     'bark': ['fibrous', 'smooth', 'mossy', 'mossy'],
     'type': ['evergreen', 'deciduous', 'deciduous', 'deciduous']}

which we've already seen.

Accessing and setting values with keys
--------------------------------------
Consistent with the behavior shown above, when accessing category values in
aggregate with keys, what is returned is a list of values for each member, in
member order::

    >>> b.categories['age']
    ['adult', 'young', 'young', 'old']

And if we access a category with a key that isn't present among all members,
``None`` is given for those members in which it's missing::

    >>> b.categories['health']
    [None, None, 'good', 'poor']

If we're interested in the values corresponding to a number of keys, we
can access these all at once with either a list::

    >>> b.categories[['health', 'bark']]
    [[None, None, 'good', 'poor'], ['fibrous', 'smooth', 'mossy', 'mossy']]

which will give a list with the values for each given key, in order by key. Or
with a set::

    >>> b.categories[{'health', 'bark'}]
    {'bark': ['fibrous', 'smooth', 'mossy', 'mossy'],
     'health': [None, None, 'good', 'poor']}

which will give a dictionary, with keys as keys and values as values.

We can also set category values for all members as if we were working
with a single member::

    >>> b.categories['height'] = 'tall'
    >>> b.categories['height']
    ['tall', 'tall', 'tall', 'tall']

or we could set the value for each member::

    >>> b.categories['height'] = ['really tall', 'middling', 'meh', 'tall']
    >>> b.categories['height']
    ['really tall', 'middling', 'meh', 'tall']

Grouping by value
-----------------
Since for a given key a Bundle may have members with a variety of values,
it can be useful to get subsets of the Bundle as a function of value for a
given key. We can do this using the
:meth:`~datreant.core.metadata.AggCategories.groupby` method::

    >>> b.categories.groupby('type')
    {'deciduous': <Bundle(['maple', 'oak', 'elm'])>,
     'evergreen': <Bundle(['sequoia'])>}

In grouping by the 'type' key, we get back a dictionary with the values present
for this key as keys and Bundles giving the corresponding members as values. We
could iterate through this dictionary and apply different operations to each
Bundle based on the value. Or we could extract out only the subset we want,
perhaps just the 'deciduous' Treants::

    >>> b.categories.groupby('type')['deciduous']
    <Bundle(['maple', 'oak', 'elm'])>

We can also group by more than one key at once::

    >>> b.categories.groupby(['type', 'health'])
    {('good', 'deciduous'): <Bundle(['oak'])>,
     ('poor', 'deciduous'): <Bundle(['elm'])>}

Now the keys of the resulting dictionary are tuples of value combinations for
which there are members. The resulting Bundles don't include some members since
not every member has both the keys 'type' and 'health'.

See the API reference for :meth:`~datreant.core.metadata.AggCategories.groupby`
for more details on how to use this method.


Operating on members in parallel
================================
Although it's common to iterate through the members of a Bundle to perform
operations on them individually, this approach can often be put in terms
of mapping a function to each member independently. A Bundle has a ``map``
method for exactly this purpose::

    >>> b.map(lambda x: (x.name, set(x.tags)))
    [('sequoia', {'huge', 'plant'}),
     ('maple', {'furniture', 'plant', 'syrup'}),
     ('oak', {'building', 'for building', 'plant'}),
     ('elm', {'building', 'firewood', 'paper', 'plant', 'shady'})]

This example isn't the most useful, but the point is that we can apply any
function across all members without much fanfare, with the results returned
in a list and in member order.

The :meth:`~datreant.core.Bundle.map` method also features a ``processes``
parameter, and setting this to an integer greater than 1 will use the
:mod:`multiprocessing` module internally to map the function across all members
using multiple processes.  For this to work, we have to give our function an
actual name so it can be serialized (pickled) by :mod:`multiprocessing`::

    >>> def get_tags(treant):
    ...     return (treant.name, set(treant.tags))
    >>> b.map(get_tags, processes=2)
    [('sequoia', {'huge', 'plant'}),
     ('maple', {'furniture', 'plant', 'syrup'}),
     ('oak', {'building', 'for building', 'plant'}),
     ('elm', {'building', 'firewood', 'paper', 'plant', 'shady'})]

For such a simple function and only four Treants in our Bundle, it's unlikely
that the parallelism gave any advantage here. But functions that need to do
more complicated work with each Treant and the data stored within its tree can 
gain much from process parallelism when applied to a Bundle of many Treants.

See the API reference for :meth:`~datreant.core.Bundle.map` for more details on
how to use this method.

API Reference: Bundle
=====================
See the :ref:`Bundle_api` API reference for more details.
