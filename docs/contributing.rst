
.. _Contributing:

========================
Contributing to datreant
========================
**datreant** is an open-source project, with its development driven by the needs
of its users. Anyone is welcome to contribute to the project, which centers
around the `datreant <https://github.com/datreant>`_ GitHub organization.

Development model
=================
datreant follows the `development model outlined by Vincent Driessen
<http://nvie.com/posts/a-successful-git-branching-model/>`_ , with the
``develop`` branch being the unstable focal point for development. The
``master`` branch merges from the develop branch only when all tests are
passing, and usually only before a release. In general, master should be usable
at all times, while develop may be broken at any particular moment.

.. _Development_env:

Setting up your development environment
=======================================
We recommend using `virtual environments
<https://pypi.python.org/pypi/virtualenv>`_ with `virtualenvwrapper
<http://virtualenvwrapper.readthedocs.org/en/latest/>`_. 
First, clone the repository::

    git clone git@github.com:datreant/datreant.git

Make a new virtualenv called ``datreant`` with::

    mkvirtualenv datreant

and make a development installation of ``datreant`` with::

    cd datreant
    pip install -e .

The ``-e`` flag will cause pip to call setup with the ``develop`` option. This
means that any changes on the source code will immediately be reflected in your
virtual environment. 

Running the tests locally
=========================
As you work on datreant, it's important to see how your changes
affected its expected behavior. With your virtualenv enabled::

    workon datreant

switch to the top-level directory of the package and run::

    py.test --cov src/ --pep8 src/
    
This will run all the tests (inside ``src/datreant/tests``), with
`coverage <https://pypi.python.org/pypi/pytest-cov>`_ and `PEP8
<https://pypi.python.org/pypi/pytest-pep8>`_ checks.

Note that to run the tests you will need to install ``py.test`` and the
coverage and PEP8 plugins into your virtualenv::

    pip install pytest pytest-cov pytest-pep8
