#! /usr/bin/python
"""Setuptools-based setup script for datreant.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(name='datreant.core',
      version='0.6.0-dev',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      packages=find_packages('src'),
      namespace_packages=['datreant'],
      package_dir={'': 'src'},
      scripts=[],
      license='BSD',
      long_description=open('README.rst').read(),
      install_requires=['asciitree', 'pathlib', 'scandir', 'six']
      )
