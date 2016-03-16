.. _tags_categories:

==========================
Differentiating Treants
==========================
Treants can be used to develop "fire-and-forget" analysis routines. Large
numbers of Treants can be fed to an analysis routine, with individual Treants
handled according to their characteristics. To make it possible to write code
that tailors its approach according to the Treant it encounters, we can use
tags and categories.

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

Categories are key-value pairs. They are particularly useful as switches for
analysis code. For example, if we have Treants with different shades of bark
(say, "dark" and "light"), we can make a category that reflects this. In this
case, we categorize ``sprout`` as "dark" ::
    
    >>> s.categories['roots'] = 'dark'
    >>> s.categories
    <Categories({'roots': 'shallow'})>

Perhaps we've written some analysis code that will take both "dark" and "light"
Treants as input but needs to handle them differently. It can see what variety
of **Treant** it is working with using ::

    >>> s.categories['roots']
    'shallow'

The keys for categories must be strings, but the values may be strings, numbers
(floats, ints), or booleans (``True``, ``False``). ``None`` may not be used since
this is used in aggregations (see :ref:`Bundles`) to indicate keys that are
absent.

Reference: Tags
===============
The class :class:`datreant.core.limbs.Tags` is the interface used by Treants to
access their tags. 

.. autoclass:: datreant.core.limbs.Tags
    :members:
    :inherited-members:

Reference: Categories
=====================
The class :class:`datreant.core.limbs.Categories` is the interface used by
Treants to access their categories.

.. autoclass:: datreant.core.limbs.Categories
    :members:
    :inherited-members:
