==========================
Differentiating Containers
==========================
**Sims** and **Groups** can be used to develop "fire-and-forget" analysis
routines. Large numbers of Containers can be fed to an analysis code to give
that code access to all trajectory and intermediate data, with individual
Containers handled according to their characteristics. To make it possible to
write code that tailors its approach according to the Container it encounters,
we can use tags and categories.

Tags are individual strings that describe a Container. Using our **Sim**
``marklar`` as an example, we can add many tags at once ::

    >>> from mdsynthesis import Sim
    >>> s = Sim('marklar')
    >>> s.tags.add('TIP4P', 'ADK', 'kinases', 'globular', 'equilibrium')
    >>> s.tags
    <Tags(['ADK', 'TIP4P', 'equilibrium', 'globular', 'kinases'])>

They can be iterated through as well ::

    >>> for tag in s.tags:
    >>>     print tag
    kinases
    globular
    ADK
    TIP4P
    equilibrium

Categories are key-value pairs of strings. They are particularly useful as
switches for analysis code. For example, if we are simulating two different
states of a protein (say, "open" and "closed"), we can make a category that
reflects this. In this case, we categorize ``marklar`` as "open" ::
    
    >>> s.categories['state'] = 'open'
    >>> s.categories
    <Categories({'state': 'open'})>

Perhaps we've written some analysis code that will take both "open" and "closed"
simulation trajectories as input but needs to handle them differently. It can
see what variety of **Sim** it is working with using ::

    >>> s.categories['state']
    'open'

Future: Querying
================
Tags and categories are two elements of Containers that will be :doc:`queryable
<Coordinator>`.

Reference: Tags
===============
The class :class:`mdsynthesis.core.aggregators.Tags` is the interface used
by Containers to access their tags. It is not intended to be used on its own,
but is shown here to give a detailed view of its methods.

.. autoclass:: mdsynthesis.core.aggregators.Tags
    :members:
    :inherited-members:

Reference: Categories
=====================
The class :class:`mdsynthesis.core.aggregators.Categories` is the interface
used by Containers to access their categories. It is not intended to be used on
its own, but is shown here to give a detailed view of its methods.

.. autoclass:: mdsynthesis.core.aggregators.Categories
    :members:
    :inherited-members:
