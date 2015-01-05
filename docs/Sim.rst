==================================
Using Sims to dissect trajectories
==================================

**Sim** objects are designed to store datasets that are obtained from a single
simulation, and they give a direct interface to trajectory data by way of the
`MDAnalysis <http://mdanalysis.googlecode.com>`__ **Universe** object.

To generate a **Sim** from scratch, we need only give it a name. This will be used
to distinguish the **Sim** from others, though it need not be unique. We can
also give it a topology and/or trajectory files as we would to an MDAnalysis
**Universe** ::
    
    >>> s = Sim('scruffy', universe=['path/to/topology', 'path/to/trajectory'])

This will create a directory ``scruffy`` that contains a single file
(``Sim.h5``).  That file is a persistent representation of the **Sim** on disk.
We can access trajectory data by way of ::

    >>> s.universe
    <Universe with 47681 atoms>

The **Sim** can also store selections by giving the usual inputs to
``Universe.selectAtoms`` ::

    >>> s.selections.add('backbone', 'name CA', 'name N', 'name C')

And the **AtomGroup** can be conveniently obtained with ::

    >>> s.selections['backbone']
    <AtomGroup with 642 atoms>

.. note:: Only selection strings are stored, not the resulting atoms of those
          selections. This means that if the topology on disk is replaced
          or altered, the results of particular selections may change.

Multiple Universes
==================
Often it is necessary to post-process a simulation trajectory to get it into a
useful form for analysis. This may involve coordinate transformations that
center on a particular set of atoms or fit to a structure, removal of water,
skipping of frames, etc. This can mean that for a given simulation, multiple
versions of the raw trajectory may be needed.

For this reason, a **Sim** can store multiple **Universe** definitions. To add
a definition, we need a topology and a trajectory file ::

    >>> s.universes.add('anotherU', 'path/to/topology', 'path/to/trajectory')
    >>> s.universes
    <Universes(['anotherU', 'main'])>

and we can make this the active **Universe** with ::

    >>> s.universes['anotherU']
    >>> s
    <Sim: 'scruffy' | active universe: 'anotherU'>

Only a single **Universe** may be active at a time. Atom selections that are
stored correspond to the currently active **Universe**, since different
selection strings may be required to achieve the same selection under a
different **Universe** definition. For convenience, we can copy the selections
corresponding to another **Universe** to the active **Universe** with ::

    >>> s.selections.copy('main')

Need two **Universe** definitions to be active at the same time? Re-generate a
second **Sim** instance from its representation on disk and activate the desired
**Universe**.

Resnums can also be stored
==========================
Depending on the simulation package used, it may not be possible to have the
resids of the protein match those given in, say, the canonical PDB structure.
This can make selections by resid cumbersome at best. For this reason, residues
can also be assigned resnums.

For example, say the resids for the protein in our **Universe** range from 1 to 214,
but they should actually go from 10 to 223. If we can't change the topology to reflect
this, we could set the resnums for these residues to the canonical values ::

    >>> prot = s.universe.selectAtoms('protein')
    >>> prot.residues.set_resnum(prot.residues.resids() + 9)
    >>> prot.residues.resnums()
    array([ 10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
            23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,
            36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,
            49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,
            62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,
            75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
            88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100,
           101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
           114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
           127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
           140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
           153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165,
           166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
           179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
           192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204,
           205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217,
           218, 219, 220, 221, 222, 223])

We can now select residue 95 from the PDB structure with ::

    >>> s.universe.selectAtoms('protein and resnum 95')

and we might save selections using resnums as well. However, resnums aren't
stored in the topology, so to avoid having to reset resnums manually each time
we load the **Universe**, we can just store the resnum definition with ::

    >>> s.universes.resnums('main', s.universe.residues.resnums())

and the resnum definition will be applied to the **Universe** both now and every
time it is activated.

Reference: Sim
==============
.. autoclass:: MDSynthesis.Sim
    :members:
    :inherited-members:

Reference: Universes
====================
The class :class:`MDSynthesis.Core.Aggregators.Universes` is the interface used
by a **Sim** to manage **Universe** definitions. It is not intended to be used
on its own, but is shown here to give a detailed view of its methods.

.. autoclass:: MDSynthesis.Core.Aggregators.Universes
    :members:
    :inherited-members:

Reference: Selections
=====================
The class :class:`MDSynthesis.Core.Aggregators.Selections` is the interface
used by a **Sim** to access its stored selections. It is not intended to be
used on its own, but is shown here to give a detailed view of its methods.

.. autoclass:: MDSynthesis.Core.Aggregators.Selections
    :members:
    :inherited-members:
