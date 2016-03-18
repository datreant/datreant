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
answering scientific questions. To make things easier, ``datreant`` is a
collection of Python packages that provide a pythonic interface to the
filesystem and the data that lives within it. It solves a boring problem, so we
can focus on interesting ones.


Stay organized
==============
``datreant`` offers a layer of flexibility and sanity to the task of analyzing data
from many studies, whether they be individual simulations or data from field
work. Its core object, the **Treant**, is designed to be subclassed: the
classes in datreant are useful on their own but vanilla by design, and are
built to be easily extended into domain-specific objects.

As an example: `MDSynthesis`_, a package for storing, recalling, and aggregating
data from molecular dynamics simulations, is built on top of datreant.

.. _`MDSynthesis`: https://github.com/datreant/MDSynthesis 


The datreant namespace
======================
``datreant`` is a namespace package, which means that it's more a package of
packages. These packages are all dependent on a central, core library, called
``datreant.core``. This documentation is for that core library.

Other packages in the ``datreant`` namespace currently include:

- `datreant.data <http://datreantdata.readthedocs.org/en/latest/>`_


Getting datreant
================
See the :doc:`installation instructions <install>` for installation details.
The package itself is pure Python, and light on dependencies by design.

If you want to work on the code, either for yourself or to contribute back to
the project, clone the repository to your local machine with::

    git clone https://github.com/datreant/datreant.core.git


Dependencies
============
The dependencies of ``datreant.core`` are light, with many being pure-Python
packages themselves. 


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

    install
    Treants
    Tags-Categories
    Trees
    Views
    Bundles
    Groups
    API


For Developers
--------------
.. toctree::
    :maxdepth: 1

    Limbs-TreeLimbs
    AggLimbs-AggTreeLimbs


Miscellaneous
-------------

.. toctree::
    :maxdepth: 1

    faq

--------------------------------------------------------------------------------

.. raw:: html

   </div>
