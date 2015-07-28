====================================
Leveraging Groups for aggregate data
====================================

A **Group** is a special type of Treant that can keep track of any number of
Treants it counts as members, and it can store datasets derived from these
objects. Just as a normal Treant manages data obtained from a single study, a
Group is designed to manage data obtained from a collection of Treants in
aggregate.

As with a normal Treant, to generate a Group from scratch, we need only give it
a name. We can also give any number of existing Treants to add them as 
members ::

    >>> from datreant import Group
    >>> g = Group('gruffy', members=[s1, s2, s3, g4, g5])
    >>> g
    <Group: 'gruffy' | 5 Members: 3 Treant, 2 Group>

This will create a directory ``gruffy`` that contains a single file
(``Group.<uuid>.h5``). That file is a persistent representation of the Group on
disk. We can access its members with ::

    >>> g.members
    <Members(['marklar', 'scruffy', 'fluffy', 'buffy', 'gorp'])>
    >>> g.members[2]
    <Treant: 'fluffy'>

and we can slice, too ::

    >>> g.members[2:]
    [<Treant: 'fluffy'>, <Group: 'buffy'>, <Group: 'gorp'>]

.. note:: Members are generated from their state files on disk upon access.
          This means that for a Group with hundreds of members, there will
          be a delay when trying to access them all at once. Once generated
          and cached, however, member access will be fast.

A Group can even be a member of itself ::

    >>> g.members.add(g)
    >>> g
    <Group: 'gruffy' | 6 Members: 3 Treant, 3 Group>
    >>> g.members[-1]
    <Group: 'gruffy' | 6 Members: 3 Treant, 3 Group>
    >>> g.members[-1].members[-1]
    <Group: 'gruffy' | 6 Members: 3 Treant, 3 Group>

As a technical aside, note that a Group returned as a member of itself
is not the same object in memory as the Group that returned it. They are
two different instances of the same Group ::

    >>> g2 = g.members[-1]
    >>> g2 is g
    False

But since they pull their state from the same file on disk, they will reflect
the same stored information at all times ::
    
    >>> g.tags.add('oaks')
    >>> g2.tags
    <Tags(['oaks'])>

Reference: Group
================
.. autoclass:: datreant.Group
    :members:
    :inherited-members:

Reference: Members
==================
The class :class:`datreant.aggregators.Members` is the interface used
by a Group to manage its members. It is not intended to be used on its own,
but is shown here to give a detailed view of its methods.

.. autoclass:: datreant.aggregators.Members
    :members:
    :inherited-members:
