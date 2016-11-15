
.. _Contributing:

========================
Contributing to datreant
========================
**datreant** is an open-source project, with its development driven by the needs
of its users. Anyone is welcome to contribute to any of its packages, which
can all be found under the `datreant <https://github.com/datreant>`_ GitHub
organization.

Development model
=================
datreant subpackages follow the `development model outlined by Vincent
Driessen <http://nvie.com/posts/a-successful-git-branching-model/>`_ , with the
``develop`` branch being the unstable focal point for development. The
``master`` branch merges from the develop branch only when all tests are
passing, and usually only before a release. In general, master should be usable
at all times, while develop may be broken at any particular moment.

.. _Development_env:

Setting up your development environment
=======================================
We recommend using `virtual environments
<https://pypi.python.org/pypi/virtualenv>`_ with `virtualenvwrapper
<http://virtualenvwrapper.readthedocs.org/en/latest/>`_. Since
datreant is a collection of subpackages, you will need to clone
whichever repositories you are interested in contributing to. Since all
datreant subpackages depend on `datreant.core
<https://github.com/datreant/datreant.core>`_, you will probably want to clone
this one::

    git clone git@github.com:datreant/datreant.core.git

Make a new virtualenv called ``datreant`` with::

    mkvirtualenv datreant

and make a development installation of ``datreant.core`` with::

    cd datreant.core
    pip install -e .

The ``-e`` flag will cause pip to call setup with the ``develop`` option. This
means that any changes on the source code will immediately be reflected in your
virtual environment. 

Repeat the same operations for any other datreant subpackage you wish to work
on. Note that other subpackages are likely to have heavier dependencies.

Running the tests locally
=========================
As you work on a datreant subpackage, it's important to see how your changes
affected its expected behavior. With your virtualenv enabled::

    workon datreant

switch to the top-level directory of the subpackage you are working on and
run::

    py.test --cov src/ --pep8 src/
    
This will run all the tests for that subpackage (often inside
``src/datreant/<subpackage_name>/tests``), with `coverage
<https://pypi.python.org/pypi/pytest-cov>`_ and `PEP8
<https://pypi.python.org/pypi/pytest-pep8>`_ checks.

Note that to run the tests you will need to install ``py.test`` and the
coverage and PEP8 plugins into your virtualenv::

    pip install pytest pytest-cov pytest-pep8
