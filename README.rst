==========================================================================
MDSynthesis: a persistence engine for intermediate molecular dynamics data
==========================================================================

|docs|

Although the raw data for any study involving molecular dynamics simulations are
the full trajectories themselves, often we are most interested in
lower-dimensional measures of what is happening. These measures may be as simple
as the distance between two specific atoms, or as complex as the percentage of
contacts relative to some native structure. In any case, it may be time-consuming
to obtain these lower-dimensional intermediate data, and so it is useful to store
them.

Stay organized
==============
MDSynthesis is designed to perform the logistics of medium-to-large-scale
analysis of many trajectories, individually or as entire groups. It should
allow the scientist to operate at a high level when working with the data,
while MDSynthesis handles the details of storing and recalling this data.

Efficiently store intermediate data for individual simulations for easy recall
------------------------------------------------------------------------------
For a given simulation trajectory, MDSynthesis gives an interface (the **Sim**
object) to the simulation data itself through `MDAnalysis`_. Data structures
generated from raw trajectories (pandas objects, numpy arrays, or any pure
python structure) can then be stored and easily recalled later. Under the hood,
datasets are stored in the efficient HDF5 format when possible.

.. _MDAnalysis: http://mdanalysis.googlecode.com

Collect aggregated data and keep track of it, too
-------------------------------------------------
**Sim** objects can be gathered into arbitrary collections with **Group** objects.
Groups can store datasets obtained from these collections, and can even contain
other Groups as members.

Query for simulation results instead of manually hunting for them
-----------------------------------------------------------------
**Note**: This feature is planned, but not yet present in the codebase.

**Sim** and **Group** objects persistently store their data to disk automatically,
but it can be tedious to navigate around the filesystem to recall them later.
The **Coordinator** object gives a single interface for querying all **Sim**
and **Group** objects it is made aware of, allowing retrieval of specific
datasets with a single line of code.

Documentation
=============
A brief user guide is available on `Read the Docs
<http://mdsynthesis.readthedocs.org/>`__.

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

.. |docs| image:: https://readthedocs.org/projects/mdsynthesis/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://readthedocs.org/projects/mdsynthesis

.. |build| image:: https://travis-ci.org/dotsdl/MDSynthesis.svg?branch=develop
    :alt: Build Status
    :scale: 100%
    :target: https://travis-ci.org/dotsdl/MDSynthesis
