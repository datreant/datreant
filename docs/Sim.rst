===========================
Using Sims to organize data
===========================

A Sim is a Container with all the machinery required to handle trajectories and
the data generated from them in an organized fashion.

To generate a Sim from scratch, we need only give it a name. This will be used
to distinguish the Sim from other Sims, though it need not be unique. We can
also give it a topology and/or trajectory files as we would to an MDAnalysis
Universe ::
    
    s = Sim('fluffy', universe=[topology, trajectory])

This will create a directory ``name`` that contains a single file (``Sim.h5``).
That file is a persistent representation of the Sim on disk. We can access
trajectory data by way of an MDAnalysis Universe ::

    s.universe

It can also store selections by giving the usual inputs to
``Universe.selectAtoms`` ::

    s.selections.add('backbone', ['name CA', 'name C', 'name O1', 'name O2'])

And the AtomGroup can be conveniently obtained with ::

    s.selections['backbone']

The Sim can also store custom data structures. These can be pandas objects
(e.g. Series, DataFrame, Panel), numpy arrays, or other python objects ::

    a = np.random.randn(100, 100)
    s.data.add('randomdata', a)

This can be recalled later with ::

    s.data['randomdata']

The real strength of the Sim is how it stores its information. Generating an
object from scratch stores the information needed to re-generate it in the
filesystem. To generate another instance of the same Sim, simply give the
directory where the state file lives ::

    s2 = Sim('fluffy/')

This Sim instance will give access to the universe, stored selections, and
stored data as before.

Reference
=========
.. autoclass:: MDSynthesis.Sim
    :members:
    :inherited-members:
