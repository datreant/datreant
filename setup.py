#! /usr/bin/python
"""Setuptools-based setup script for datreant.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='datreant',
      version='0.6.0-dev',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      packages=[
          'datreant',
          'datreant.backends',
          'datreant.tests'],
      scripts=[],
      license='BSD',
      long_description=open('README.rst').read(),
      install_requires=['scandir', 'six']
      )
