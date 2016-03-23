.. _tags_categories:

==========================
Differentiating Treants
==========================
Treants can be used to develop "fire-and-forget" analysis routines. Large
numbers of Treants can be fed to an analysis routine, with individual Treants
handled according to their characteristics. To make it possible to write code
that tailors its approach according to the Treant it encounters, we can use
tags and categories.


.. _Tags_guide:

Using tags
==========
Tags are individual strings that describe a Treant. Using our Treant
``sprout`` as an example, we can add many tags at once ::

    >>> from datreant import Treant
    >>> s = Treant('sprout')
    >>> s.tags.add('elm', 'mirky', 'misty')
    >>> s.tags
    <Tags(['elm', 'mirky', 'misty'])>

They can be iterated through as well ::

    >>> for tag in s.tags:
    >>>     print tag
    elm
    mirky
    misty

Or checked for membership ::

    >>> 'mirky' in s.tags
    True

Tags function as sets
---------------------
Since the tags of a Treant behave as a set, we can do set operations directly,
such as subset comparisons::

    >>> {'elm', 'misty'} < s.tags
    True

unions::

    >>> {'moldy', 'misty'} | s.tags
    {'elm', 'mirky', 'misty', 'moldy'}

intersections::

    >>> {'elm', 'moldy'} & s.tags
    {'elm'}

differences::

    >>> s.tags - {'moldy', 'misty'}
    {'elm', 'mirky'}

or symmetric differences::

    >>> s.tags ^ {'moldy', 'misty'}
    {u'elm', u'mirky', 'moldy'}

It is also possible to set the tags directly::

    >>> s.tags = s.tags | {'moldy', 'misty'}
    >>> s.tags
    <Tags(['elm', 'mirky', 'misty', 'moldy'])>

API Reference: Tags
-------------------
See the :ref:`Tags_api` API reference for more details.


.. _Categories_guide:

Using categories
================
Categories are key-value pairs. They are particularly useful as switches for
analysis code. For example, if we have Treants with different shades of bark
(say, "dark" and "light"), we can make a category that reflects this. In this
case, we categorize ``sprout`` as "dark" ::
    
    >>> s.categories['bark'] = 'dark'
    >>> s.categories
    <Categories({'bark': 'dark'})>

Perhaps we've written some analysis code that will take both "dark" and "light"
Treants as input but needs to handle them differently. It can see what variety
of **Treant** it is working with using ::

    >>> s.categories['bark']
    'dark'

The keys for categories must be strings, but the values may be strings, numbers
(floats, ints), or booleans (``True``, ``False``). 

.. note:: ``None`` may not be used as a category value since this is used in
          aggregations (see :ref:`Bundles`) to indicate keys that are absent.

API Reference: Categories
-------------------------
See the :ref:`Categories_api` API reference for more details.


Filtering and grouping on tags and categories
=============================================
Tags and categories are especially useful for filtering and grouping Treants.
See :ref:`Bundles` for the details on how to flexibly do this.
