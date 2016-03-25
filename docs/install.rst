========================
Installing datreant.core
========================
You can install ``datreant.core`` from `PyPI <https://pypi.python.org/>`_ using pip::

    pip install datreant.core

It is also possible to use ``--user`` to install into your user's site-packages
directory::

    pip install --user datreant.core

All datreant packages currently support the following Python versions::

- 2.7
- 3.3
- 3.4
- 3.5


Dependencies
============
The dependencies of ``datreant.core`` are light, with many being pure-Python
packages themselves. The current dependencies are::

- asciitree
- pathlib
- scandir
- six
- fuzzywuzzy

These are automatically installed when installing ``datreant.core``.

Installing from source
======================

To install from source, clone the repository and switch to the master branch ::

    git clone git@github.com:datreant/datreant.core.git
    cd datreant.core
    git checkout master

Installation of the packages is as simple as ::

    pip install .

This installs ``datreant.core`` in the system wide python directory; this may
require administrative privileges. If you have a virtualenv active, it will
install the package within your virtualenv. See :ref:`Development_env` for more
on setting up a proper development environment.

It is also possible to use ``--user`` to install into your user's site-packages
directory::

    pip install --user .
