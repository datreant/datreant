#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='mdsynthesis',
      version='0.4.0-dev',
      author='David Dotson', 
      author_email='dotsdl@gmail.com',
      packages=['mdsynthesis', 'mdsynthesis.core'],
      scripts=['scripts/mds_convert_pre_to_0.5.0.py'],
      license='GPL 2',
      long_description=open('README.rst').read(),
      requires=['pandas', 'tables', 'h5py', 'MDAnalysis', 'scandir']
     )
