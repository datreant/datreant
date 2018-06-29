#! /usr/bin/python
"""Setuptools-based setup script for datreant.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(
    name='datreant',
    version='1.0.0',
    description='persistent, pythonic trees for heterogeneous data',
    author='David Dotson',
    author_email='dotsdl@gmail.com',
    url='http://datreant.org/',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    scripts=[],
    entry_points={'console_scripts':
                  ['datreant_07to1=datreant.scripts.datreant_07to1.main']},
    license='BSD',
    long_description=open('README.rst').read(),
    tests_require=['numpy', 'pytest>=2.10', 'mock'],
    install_requires=[
        'asciitree', 'pathlib2', 'scandir', 'six', 'fuzzywuzzy',
        'python-Levenshtein', 'pyparsing'])
