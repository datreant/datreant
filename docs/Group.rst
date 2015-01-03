====================================
Leveraging Groups for aggregate data
====================================

A Group is a Container that keeps track of any number of Sims and Groups added
to it as members, and it can store datasets derived from these objects in the
same way as Sims.

To generate a Group object from scratch, we need only give it a name. We
can also give any number of Sim and/or Groups as an argument ::

    g = Group('gruffy', members=[s1, s2, s3, g4, g5])

The Group can store custom data structures the same way as Sims. These can
be pandas objects (e.g. Series, DataFrame, Panel), numpy arrays, or other
python objects ::

    a = np.random.randn(100, 100)
    g.data.add('randomdata', a)

This can be recalled later with ::

    g.data['randomdata']

The real strength of the Group is how it stores its information. Generating an
object from scratch stores the information needed to re-generate it in the
filesystem. To generate another instance of the same Group, simply give the
directory where the state file lives ::

    g2 = Group('gruffy/')

This Group instance will give access to its members and stored data as before.

Reference
=========
.. autoclass:: MDSynthesis.Group
    :members:
    :inherited-members:
