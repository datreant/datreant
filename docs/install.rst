============
Installation
============
There are no official releases of MDSynthesis yet, but the development
branch on GitHub gives the most current state of the package. Clone
the repository and then switch to this branch ::

    git clone git@github.com:dotsdl/MDSynthesis.git
    cd MDSynthesis
    git checkout develop

Installation is as simple as ::

    python setup.py build
    python setup.py install

This installs MDSynthesis in the system wide python directory; this may
require administrative privileges.

It is also possible to use ``--prefix``, ``--home``, or ``--user`` options for
setup.py to install in a different (probably your private) python directory
hierarchy. ``python setup.py install --help`` should show you your options.

