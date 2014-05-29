# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDSynthesis 

"""
:mod:`MDSynthesis` --- an organizational framework using MDAnalysis
====================================================================

MDSynthesis is designed to address the logistical aspects of molecular dynamics
trajectory analysis. Whereas MDAnalysis gives the computational tools to dissect
trajectories, MDSynthesis provides a framework for automatically organizing the
results. This allows you (the scientist) to focus on your science, letting the
computer handle the lower-level logistical details. 

Generating Sim objects
----------------------

In an interactive python session, import the package::

  >>> import MDSynthesis

The Sim class is immediately available. You can create a Sim object by giving it
a MDAnalysis Universe::

  >>> u = MDAnalysis.Universe(PSF, DCD)
  >>> s = MDSynthesis.Sim(u, location=".")

The Sim object is now loaded in the python session, but it is also stored in
the filesystem. In this case, we specified that it be stored in the current
directory. You can find this object in `./MDSynthesis/Sim/`, where it is
really just a single YaML file that stores its metadata.

Now, leave the session. Starting a new session, you can reload the object by
simply specifiying the directory in which the metadata file is stored::

  >>> import MDSynthesis
  >>> s = MDSynthesis.Sim('./MDSynthesis/Sim/')

The metadata is loaded, and the Sim object is identical to the one initialized
before using a universe.

Using Analysis objects
----------------------

The real strength of MDSynthesis is in its analysis framework. The Analysis
class in :mod:`MDSynthesis.Operators` can be used as a base class for any
timeseries analysis. If we built an analysis object called `Salt` that looked
at ion interactions with a protein, we could then run this on our simulation object
with::

  >>> salt = Salt(s) 
  >>> salt.run()

For this Sim object, the resulting data is pickled and stored in 
`./MDSynthesis/Sim/Salt/Salt.pkl`. It can be reloaded for immediate parsing
with::

  >>> s.load('Salt')

and unloaded with::

  >>> s.unload('Salt')

This allows one to collect large numbers of unrelated datasets for a specific
simulation system, then easily recall select datasets for analysis later. All
the while, the datasets remain organized.

Further automation
------------------

Instead of specifying a location to store the Sim object as before, we can
specify the project directory out of which we are working::

  >>> u = MDAnalysis.Universe(PSF, DCD)
  >>> s = MDSynthesis.Sim(u, projectdir="/home/monty/Projects/Spam/Eggs/")

If the trajectory file is stored in `/home/monty/Projects/Spam/Eggs/WORK/system1`, 
then the Sim object will be placed in 
`/home/monty/Projects/Spam/Eggs/MDSynthesis/Sim/WORK/system1`. One can remove
portions of the directory structure after the object name (`Sim`) using the
`pluck_segment` keyword.

The advantage with this is that, once created, the Sim object can be reloaded
from the directory in which it is stored. You need never track down the
structure file or trajectory file ever again. 

The Sim class can also be used as a base class for a container that is more
specific to the type of simulations one is working with. All classes
in MDSynthesis are designed with easy inheritance in mind.

.. SeeAlso:: :class:`MDSynthesis.Containers.Sim` for details


"""

__version__ = "0.3.0-dev"  # NOTE: keep in sync with RELEASE in setup.py

__all__ = ['Sim', 'Group', 'Analysis', 'MetaAnalysis']

# Bring some often used objects into the current namespace
from Database import Database
from Containers import Sim, Group, SuperGroup
from Operators import Analysis, MetaAnalysis, MacroAnalysis
