.. datreant documentation master file

===========================================================
datreant: persistent, pythonic trees for heterogeneous data
===========================================================
In many fields of science, especially those analyzing experimental or
simulation data, there is often an existing ecosystem of specialized tools and 
file formats which new tools must work around, for better or worse.
Furthermore, centralized database solutions may be suboptimal for data
storage for a number of reasons, including insufficient hardware
infrastructure, variety and heterogeneity of raw data, the need for data
portability, etc. This is particularly the case for fields centered around
simulation: simulation systems can vary widely in size, composition, rules,
paramaters, and starting conditions. And with increases in computational power,
it is often necessary to store intermediate results obtained from large amounts
of simulation data so it can be accessed and explored interactively.

These problems make data management difficult, and serve as a barrier to
answering scientific questions. To make things easier, **datreant** is a Python
package that addresses the tedious and time-consuming logistics of intermediate
data storage and retrieval. It solves a boring problem, so we can focus on
interesting ones.

Stay organized
==============
datreant is offers a layer of flexibility and sanity to the task of analyzing
data from many studies, whether they be individual simulations or data from
field work. It is a library that is designed to be subclassed: the classes in
datreant are useful on their own but vanilla by design, and are built to be
easily extended into domain-specific objects.

As an example: `MDSynthesis`_, a package for storing, recalling, and aggregating
data from molecular dynamics simulations, is built on top of datreant.

.. _`MDSynthesis`: https://github.com/datreant/MDSynthesis 

Efficiently store intermediate data from individual studies for easy recall
---------------------------------------------------------------------------
For handling data from a single study, datreant gives the **Treant** object.  A
**Treant** can store data structures generated from raw data (pandas objects,
numpy arrays, or any pure python structure) for easy recall later. Under the
hood, datasets are stored in the efficient `HDF5`_ format when possible.

.. _`HDF5`: https://www.hdfgroup.org/HDF5/whatishdf5.html

Collect aggregated data and keep track of it, too
-------------------------------------------------
**Treants** can be gathered into arbitrary collections with **Group** objects.
**Groups** can store datasets obtained from these collections, and can even
contain other **Groups** as members. **Groups** can keep track of any
**Treant**-derived subclasses, even domain-specific ones that you've defined
for getting your work done.

Query for study results instead of manually hunting for them
------------------------------------------------------------
**Note**: This feature is planned, but not yet present in the codebase.

**Treant** and **Group** objects persistently store their data to disk
automatically, but it can be tedious to navigate around the filesystem to
recall them later.  The **Coordinator** object gives a single interface for
querying all **Treants** and **Groups** it is made aware of, allowing retrieval
of specific datasets with a single line of code.

Getting datreant
================
See the :doc:`installation instructions <install>` for installation details.
The package itself is pure Python, but it is dependent on `HDF5`_ libraries
and the Python interfaces to these.

If you want to work on the code, either for yourself or to contribute back to
the project, clone the repository to your local machine with::

    git clone https://github.com/datreant/datreant.git

Dependencies
============
* `pandas`_: 0.16.1 or higher
* `PyTables`_: 3.2.0 or higher
* `h5py`_: 2.5.0 or higher
* `scandir`_: 1.0 or higher

.. _`pandas`: http://pandas.pydata.org/
.. _`PyTables`: http://www.pytables.org/
.. _`h5py`: http://www.h5py.org/
.. _`scandir`: https://pypi.python.org/pypi/scandir

Contributing
============
This project is still under heavy development, and there are certainly rough
edges and bugs. Issues and pull requests welcome! Check out our `contributor's guide`_
if you learn how to get started with contributing back.

.. _`contributor's guide`: https://github.com/datreant/datreant/wiki/Contributing

.. raw:: html

    <div style="display:none">

--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
    :maxdepth: 1

    install
    Treant
    Data
    Group
    Tags-Categories
    Coordinator

Misc
----

.. toctree::
    :maxdepth: 1

    faq

--------------------------------------------------------------------------------

.. raw:: html

   </div>
