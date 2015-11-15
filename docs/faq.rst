==========================
Frequently Asked Questions
==========================

1. What are some benefits of datreant's approach to data management?

    *It's daemonless.* 
        There is no server process required on the machine to
        read/write persistent objects. This is valuable when working with data on
        remote resources, where it might not be possible to set up one's own daemon
        to talk to.

    *treants are portable.*
        Because they store their state in HDF5 (which itself is portable across
        systems, handling endianess, etc.), and because treants store their
        data in the filesystem, they are easy to move around piecemeal. If you
        want to use a treant on a remote system, but don't want to drag all
        its stored datasets with it, you can copy only what you need.

        Contrast this with many database solutions, in which you either copy the
        whole database somehow, or slurp the pieces of data out that you want.
        Most database solutions can be rather slow to do this, to my knowledge.

    *treants are independent.*
        Although Groups are aware of their members, treants work
        independently from one another. If you want to use only basic
        treants or Sims, that works just fine. If you want to use Groups,
        that works, too. If you want to use the Coordinator (not yet
        implemented; think a thin database that treants share their info
        with so they can be quickly queried from one place in the filesystem),
        then you can, but you don't have to. You don't have to buy the whole
        farm to ride the horse, in a sense.

    *treants have a structure in the filesystem.* 
        This means that all the shell tools we know and love are available to
        work with their contents, which might include plaintext files, figures,
        topology files, trajectories, random pickles, ipython notebooks, html
        files, etc. Basically, treants are as versatile as the filesystem
        is, at least when it comes to storage.

2. What are some disadvantages of datreant's design?
    
    *treants could be anywhere in the filesystem.*
        This is mostly a problem for Groups, which allow aggregation of other
        treants. If a member is moved, the Group has no way of knowing where
        it went; we've built machinery to help it find its members, but these
        will always be limited to filesystem search methods (some quite good,
        but still). If these objects lived in a single database, this wouldn't
        be an issue.

    *Queries on object metadata will be slower than a central database.*
        We want Groups and Bundles (in-memory Groups, basically) to be able to
        run queries against their members' characteristics, returning subsets
        matching the query.  Since these queries have to be applied against
        these objects and not against a single table somewhere, it will be
        relatively slow. 

        The Coordinator is an answer to this problem, albeit an imperfect one.
        The idea is that you can make a Coordinator, which is a small
        daemonless database (perhaps SQLite, but could be HDF5), and you can
        add treants for it to be aware of with something like::

            import mdsynthesis as mds

            co = mds.Coordinator('camelot')
            co.add('moe', 'larry', 'curly')

            # could also let the coordinator do a downward search and add all
            # treants it finds
            co.discover()
            
        This awareness is bi-directional: a Coordinator is aware of its
        members, and its members are aware of the Coordinator, and where it
        lives. The Coordinator will store tables of member attributes for fast
        queries, and these tables will be updated by members as they themselves
        are updated. So whenever we have::

            python

            c = mds.treant('moe')
            c.categories['bowlcut'] = True
   
        the treant updates both its state file and the Coordinator(s) it is
        affiliated with. This is in contrast to Groups, of which members are
        unaware. This is by design: the idea is that treants are likely to be
        members of many different Groups all over the filesystem; there would be
        comparatively fewer Coordinators in use, which have a performance hit to a
        treant for each affiliation.

        Obvious problem: there are probably a lot of ways for Coordinators and their
        members to get out-of-sync. A single database with everything inside avoids
        this entirely.

    *File locking is less efficient for multiple-read/write under load than a smart daemon process/scheduler.* 
        The assumption we make is that treants are primarily read, and only
        occasionally written to. This is assumed for their data and their
        metadata. They are not designed to scale well if the same parts are
        being written to and read at the same time by many processes.

        Having treants exist as separate files (state files and data all
        separate) does mitigate this potential for gridlock, which is one
        reason we favor many files over few. But it's still something to be
        aware of.
