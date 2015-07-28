==========================
Differentiating Treants
==========================
Treants can be used to develop "fire-and-forget" analysis routines. Large
numbers of Treants can be fed to an analysis code to give that code access to
all stored datasets, with individual Treants handled according to their
characteristics. To make it possible to write code that tailors its approach
according to the Treant it encounters, we can use tags and categories.

Tags are individual strings that describe a Treant. Using our Treant
``sprout`` as an example, we can add many tags at once ::

    >>> from datreant import Treant
    >>> t = Treant('sprout')
    >>> t.tags.add('elm', 'mirky', 'misty')
    >>> t.tags
    <Tags(['elm', 'mirky', 'misty'])>

They can be iterated through as well ::

    >>> for tag in t.tags:
    >>>     print tag
    elm
    mirky
    misty

Categories are key-value pairs of strings. They are particularly useful as
switches for analysis code. For example, if we have Treants with different
depths of roots (say, "deep" and "shallow"), we can make a category that
reflects this. In this case, we categorize ``sprout`` as "shallow" ::
    
    >>> t.categories['roots'] = 'shallow'
    >>> t.categories
    <Categories({'roots': 'shallow'})>

Perhaps we've written some analysis code that will take both "deep" and "shallow"
Treants as input but needs to handle them differently. It can see what variety
of **Treant** it is working with using ::

    >>> t.categories['roots']
    'shallow'

Future: Querying
================
Tags and categories are two elements of Treants that will be :doc:`queryable
<Coordinator>`.

Reference: Tags
===============
The class :class:`datreant.aggregators.Tags` is the interface used
by Treants to access their tags. It is not intended to be used on its own,
but is shown here to give a detailed view of its methods.

.. autoclass:: datreant.aggregators.Tags
    :members:
    :inherited-members:

Reference: Categories
=====================
The class :class:`datreant.aggregators.Categories` is the interface
used by Treants to access their categories. It is not intended to be used on
its own, but is shown here to give a detailed view of its methods.

.. autoclass:: datreant.aggregators.Categories
    :members:
    :inherited-members:
