.. MDSynthesis documentation master file, created by
   sphinx-quickstart2 on Mon Dec 29 20:48:44 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=============================================================
MDSynthesis: a persistence engine for molecular dynamics data
=============================================================
Although the raw data for any study involving molecular dynamics simulations are
the full trajectories themselves, often we are most interested in
lower-dimensional measures of what is happening. These measures may be as simple
as the distance between two specific atoms, or as complex as the percentage of
contacts relative to some native structure. Some measures may even be
comparisons of one or more trajectories against each other. In any case, it may
be time-consuming to obtain these lower-dimensional intermediate data, and so
it is useful to store them.

Stay organized
==============
MDSynthesis is designed to perform the logistics of medium-to-large-scale
analysis of many trajectories, individually or as entire groups. It is intended
to allow the scientist to operate at a high level when working with the data,
while letting MDSynthesis handle the details of storing and recalling this
data. 

In other words, MDSynthesis lets the computer do the boring work of keeping
track of where things are and how they are stored.

Efficiently store intermediate data for individual simulations for easy recall
------------------------------------------------------------------------------
For a given simulation trajectory, MDSynthesis gives an interface (the :doc:`Sim <Sim>`
object) to the simulation data itself through `MDAnalysis`_. Data structures
generated from raw trajectories (pandas objects, numpy arrays, or any pure
python structure) can then be stored and easily recalled later. Under the hood,
datasets are stored in the efficient HDF5 format when possible.

.. _MDAnalysis: http://mdanalysis.googlecode.com

Collect aggregated data and keep track of it, too
-------------------------------------------------
:doc:`Sim <Sim>` objects can be gathered into arbitrary collections with
:doc:`Group <Group>` objects.  Groups can store datasets obtained from these
collections, and can even contain other Groups as members.

Query for simulation results instead of manually hunting for them
-----------------------------------------------------------------
.. note:: This feature is planned, but not yet present in the codebase.

:doc:`Sim <Sim>` and :doc:`Group <Group>` objects persistently store their data
to disk automatically, but it can be tedious to navigate around the filesystem
to recall them later.  The :doc:`Coordinator <Coordinator>` object gives a
single interface for querying all :doc:`Sim <Sim>` and :doc:`Group <Group>`
objects it is made aware of, allowing retrieval of specific datasets with a
single line of code.

Getting MDSynthesis
===================
We have yet to make an official release, but you can get the current state
of the codebase from the `development branch on GitHub 
<https://github.com/dotsdl/MDSynthesis>`__.

See the :doc:`installation instructions <install>` to set it up.

Dependencies
============
* MDAnalysis: 0.8.0 or higher
* Pandas: 0.15.0 or higher
* PyTables: 3.1.0 or higher
* h5py: 2.3.0 or higher

Contributing
============
This project is still under heavy development, and there are certainly rough
edges and bugs. Issues and pull requests welcome!

.. raw:: html

    <div style="display:none">

--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
    :maxdepth: 1

    data
    Sim
    Group
    tags-categories
    Coordinator
    interactive

Misc
----

.. toctree::
    :maxdepth: 1

    install
    whatsnew
    faq

--------------------------------------------------------------------------------

.. raw:: html

   </div>

