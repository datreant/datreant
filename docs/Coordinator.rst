Query and high-level control with Coordinators
==============================================
Because Treants store their information neatly in their state files, this data
can be aggregated and queried. This allows whole selections of Treants to be
manipulated without needing to hunt them down in the filesystem. The
**Coordinator** object gives an interface for doing this. Treants that are
associated with a given Coordinator will report changes to their state files as
they are made, giving the Coordinator a thin copy of all Treants it is made
aware of.

This feature is not yet implemented.

.. .. autoclass:: datreant.Coordinator
    :members:
    :inherited-members:
