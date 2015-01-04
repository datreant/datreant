Getting started
===============
MDSynthesis is not an analysis code. On its own, it does not produce output
data given raw simulation data as input. Its scope is limited to the boring
but tedious task of data management and storage. It is intended to bring
value to analysis results by making them easily accessible now and later.

As such, the basic functionality of MDSynthesis is condensed into only two
objects, sometimes referred to as *Containers* in the documentation. These are
the :doc:`Sim <Sim>` and :doc:`Group <Group>` objects.

In brief, a **Sim** is designed to manage and give access to the data corresponding
to a single simulation (the raw trajectory(s), as well as analysis results); a
**Group** gives access to any number of **Sim** or other **Group** objects
it has as members (including perhaps itself), and can store analysis results
that pertain to these members collectively. Both types of Container store
their underlying data persistently to disk on the fly. The file locking needed
for each transaction is handled automatically, so more than one python session
can be working with any number of instances of the same Container at the same
time.

Containers store their data as directory structures in the file system. Generating
a new **Sim**, for example, with the following ::
    
    [ python session 1 ]
    >>> import MDSynthesis as mds
    >>> s = mds.Sim('marklar')

creates a directory called ``marklar`` in the current working directory. It contains
a single file at the moment ::

    [ shell ]
    > ls marklar
    Sim.h5

This is the state file containing all the information needed to regenerate an
identical instance of this **Sim**. In fact, we can open a separate python
session (go ahead!) and regenerate this **Sim** immediately there ::

    [ python session 2 ]
    >>> import MDSynthesis as mds
    >>> s = mds.Sim('marklar')

Making a modification to the **Sim** in one session, perhaps by adding a tag,
will be reflected in the **Sim** in the other session ::

    [ python session 1 ]
    >>> s.tags.add('TIP3P')

    [ python session 2 ]
    >>> s.tags
    Tags(['TIP3P'])

This is because both objects pull their identifying information from the same
file on disk; they store almost nothing in memory.


