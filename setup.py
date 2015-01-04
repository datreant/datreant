#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='MDSynthesis',
      version='0.4.0-dev',
      author='David Dotson', 
      author_email='dotsdl@gmail.com',
      packages=['MDSynthesis', 'MDSynthesis.Core'],
      license='GPL 2',
      long_description=open('README.rst').read(),
      requires=['MDAnalysis', 'pandas', 'tables', 'h5py']
     )
