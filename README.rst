===========================================================
datreant: persistent, pythonic trees for heterogeneous data
===========================================================

|docs| |build| |cov|

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

For more information on what **datreant** is and what it does, check out the
`official documentation`_.

.. _`official documentation`: http://datreant.readthedocs.org/

Getting datreant
================
See the `installation instructions`_ for installation details.
The package itself is pure Python.

If you want to work on the code, either for yourself or to contribute back to
the project, clone the repository to your local machine with::

    git clone https://github.com/datreant/datreant.core.git

.. _`installation instructions`: http://datreant.readthedocs.org/en/develop/install.html

Contributing
============
This project is still under heavy development, and there are certainly rough
edges and bugs. Issues and pull requests welcome! 

Check out our `contributor's guide`_ to learn how to get started with
contributing back.

.. _`contributor's guide`: http://datreant.readthedocs.org/en/develop/contributing.html

.. |docs| image:: https://readthedocs.org/projects/datreant/badge/?version=develop
    :alt: Documentation Status
    :scale: 100%
    :target: http://datreant.readthedocs.org/en/develop/?badge=develop

.. |build| image:: https://travis-ci.org/datreant/datreant.core.svg?branch=develop
    :alt: Build Status
    :target: https://travis-ci.org/datreant/datreant.core

.. |cov| image:: http://codecov.io/github/datreant/datreant.core/coverage.svg?branch=develop
    :alt: Code Coverage
    :scale: 100%
    :target: http://codecov.io/github/datreant/datreant.core?branch=develop

