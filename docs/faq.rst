==========================
Frequently Asked Questions
==========================

1. What are some benefits of datreant's approach to data management?

    *It's daemonless.* 
        There is no server process required on the machine to
        read/write persistent objects. This is valuable when working with data on
        remote resources, where it might not be possible to set up one's own daemon
        to talk to.

    *Treants are portable.*
        Because they store their state in JSON, and because treants store their
        data in the filesystem, they are easy to move around piecemeal. If you
        want to use a treant on a remote system, but don't want to drag all its
        stored datasets with it, you can copy only what you need.

        Contrast this with many database solutions, in which you either copy the
        whole database somehow, or slurp the pieces of data out that you want.
        Most database solutions can be rather slow to do this.

    *Treants have a structure in the filesystem.* 
        This means that all the shell tools we know and love are available to
        work with their contents, which might include plaintext files, figures,
        topology files, simulation trajectories, random pickles, ipython
        notebooks, html files, etc. Basically, treants are as versatile as the
        filesystem is, at least when it comes to storage.

2. What are some disadvantages of datreant's design?

    *Queries on object metadata will be slower than a central database.*
        We Bundles to be able to run queries against their members'
        characteristics, returning subsets matching the query.  Since these
        queries have to be applied against these objects and not against a
        single table somewhere, it will be relatively slow. 

    *File locking is less efficient for multiple-read/write under load than a smart daemon process/scheduler.* 
        The assumption we make is that Treants are primarily read, and only
        occasionally written to. This is assumed for their data and their
        metadata. They are not designed to scale well if the same parts are
        being written to and read at the same time by many processes.

        Having Treants exist as separate files (state files and data all
        separate) does mitigate this potential for gridlock, which is one
        reason we favor many files over few. But it's still something to be
        aware of.
