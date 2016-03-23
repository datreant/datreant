==============================
Persistent Bundles with Groups
==============================
A **Group** is a Treant that can keep track of any number of Treants it counts
as members. It can be thought of as having a persistent Bundle attached to it,
with that Bundle's membership information stored inside the Group's state file.
Since Groups are themselves Treants, they are useful for managing data obtained
from a collection of Treants in aggregate.

.. note:: Since Groups are Treants, everything that applies to Treants applies
          to Groups as well. They have their own tags, categories, and
          directory trees, independent of those of their members.

As with a normal Treant, to generate a Group from scratch, we need only give it
a name ::

    >>> from datreant.core import Group
    >>> g = Group('forest')
    >>> g
    <Group: 'forest'>

We can also give any number of existing Treants to add them as 
members, including other Groups ::

    >>> g.members.add('sequoia', 'maple', 'oak', 'elm')
    >>> g
    <Group: 'forest' | 4 Members>

We can access its members directly with::

    >>> g.members
    <MemberBundle([<Treant: 'sequoia'>, <Treant: 'maple'>, <Treant: 'oak'>, <Treant: 'elm'>])>

which yields a **MemberBundle**. This object works exactly the same as a
Bundle, with the only difference that the order and contents of its membership
are stored within the Group's state file. Obtaining subsets of the
MemberBundle, such as by slicing, indexing, filtering on tags, etc., will
always yield a Bundle::

    >>> g.members[2:]
    <Bundle([<Treant: 'oak'>, <Treant: 'elm'>])>

.. note:: Members are generated from their state files on disk upon access.
          This means that for a Group with hundreds of members, there will
          be a delay when trying to access them all at once. Once generated
          and cached, however, member access will be fast.

A Group can even be a member of itself ::

    >>> g.members.add(g)
    >>> g
    <Group: 'forest' | 5 Members>
    >>> g.members[-1]
    <Group: 'forest' | 5 Members>


Groups find their members when they go missing
==============================================
For each member Treant, a Group stores the Treant's type ('Treant', 'Group'),
its uuid, and last known absolute and relative locations. Since Treants are
directories in the filesystem, they can very well be moved. What happens
when we load up an old Group whose members may no longer be where they once
were?

Upon access to its MemberBundle, a Group will try its best to track down its
members. It will first check the last known locations for each member, and
for those it has yet to find it will begin a downward search starting from its
own tree, then its parent tree, etc. It will continue this process until it
either finds all its members, hits the root of the filesystem, or its search
times out.

If the Group you are using fails to find some of its members before timing out,
you can set the maximum search time to a longer timeout, in seconds::

    >>> g.members.searchtime = 60

Otherwise, you will need to remove the missing members or find them yourself.


API Reference: Group
====================
See the :ref:`Group_api` API reference for more details.
