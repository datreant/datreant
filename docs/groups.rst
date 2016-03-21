====================================
Leveraging Groups for aggregate data
====================================

A **Group** is a special type of Treant that can keep track of any number of
Treants it counts as members. Just as a normal Treant can be used to manage
data obtained from a single study, a Group is useful for managing data obtained
from a collection of Treants in aggregate.

As with a normal Treant, to generate a Group from scratch, we need only give it
a name ::

    >>> from datreant.core import Group
    >>> g = Group('gruffy')
    >>> g
    <Group: 'gruffy'>

We can also give any number of existing Treants to add them as 
members, including other Groups ::

    >>> members=[s1, s2, s3, g4, g5])
    <Group: 'gruffy' | 5 Members>

This will create a directory ``gruffy`` that contains a single file
(``Group.<uuid>.json``). That file is a persistent representation of the Group
on disk. We can access its members with ::

    >>> g.members
    <MemberBundle(['marklar', 'scruffy', 'fluffy', 'buffy', 'gorp'])>
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

As an aside, note that a Group returned as a member of itself is not the same
object in memory as the Group that returned it. They are two different
instances of the same Group ::

    >>> g2 = g.members[-1]
    >>> g2 is g
    False

But since they pull their state from the same file on disk, they will reflect
the same stored information at all times ::
    
    >>> g.tags.add('oaks')
    >>> g2.tags
    <Tags(['oaks'])>

Groups find their members when they go missing
==============================================



API Reference: Group
====================
See the :ref:`Group_api` API reference for more details.
