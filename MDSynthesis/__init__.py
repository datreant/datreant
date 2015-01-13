# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDSynthesis

"""
:mod:`MDSynthesis` --- a persistence engine for molecular dynamics data
========================================================================

MDSynthesis is designed to address the logistical aspects of molecular dynamics
trajectory analysis. Whereas MDAnalysis gives the computational tools to
dissect trajectories, MDSynthesis provides a framework for automatically
organizing the results. This allows you (the scientist) to focus on your
science, letting the computer handle the lower-level logistical details.

.. SeeAlso:: :class:`MDSynthesis.Containers.Sim` 
             :class:`MDSynthesis.Containers.Group`


"""

__version__ = "0.4.0-dev"  # NOTE: keep in sync with RELEASE in setup.py

__all__ = ['Sim', 'Group', 'Coordinator']

# Bring some often used objects into the current namespace
#from Coordinator import Coordinator
from Containers import Sim, Group
from Coordinator import Coordinator
import Core
from Core.Workers import Bundle
