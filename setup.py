#! /usr/bin/python
"""Setuptools-based setup script for datreant.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='datreant',
      version='0.5.1',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      packages=['datreant', 'datreant.tests'],
      scripts=[],
      license='BSD',
      long_description=open('README.rst').read(),
      install_requires=[
          'numpy',
          'pandas',
          'tables',
          'h5py',
          'scandir',
          'PyYAML'
          ]
      )
