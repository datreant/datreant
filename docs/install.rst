========================
Installing datreant.core
========================
There are no official releases of datreant yet, but the master branch on GitHub
gives the most current state of the package. 

.. note:: Python 2.7 or above is required for all datreant packages.

To install from source, clone the repository and switch to the master branch ::

    git clone git@github.com:datreant/datreant.core.git
    cd datreant.core
    git checkout master

Installation of the packages is as simple as ::

    python setup.py install

This installs datreant in the system wide python directory; this may require
administrative privileges.

It is also possible to use ``--prefix``, ``--home``, or ``--user`` options for
setup.py to install in a different (probably your local user) Python
site-packages directory. ``python setup.py install --help`` should show you
your options.
