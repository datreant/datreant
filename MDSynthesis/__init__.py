# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDSynthesis 

"""
:mod:`MDSynthesis` --- an organizational framework using MDAnalysis
====================================================================


.. rubric:: Included algorithms

If you use the RMSD alignment code (:mod:`MDAnalysis.analysis.align`) that
uses the :mod:`~MDAnalysis.core.qcprot` module please also cite

    Douglas L. Theobald. Rapid calculation of RMSD using a
    quaternion-based characteristic polynomial. Acta Crystallographica
    A 61 (2005), 478-480.

    Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald. Fast
    determination of the optimal rotational matrix for macromolecular
    superpositions. J. Comput. Chem. 31 (2010), 1561-1563.

Getting started
---------------

Import the package::

  >>> import MDSynthesis

(note that not everything in MDAnalysis is imported right away; for
additional functionality you might have to import sub-modules
separately, e.g. for RMS fitting ``import MDAnalysis.analysis.align``.)

Build a "universe" from a topology (PSF, PDB) and a trajectory (DCD, XTC/TRR);
here we are assuming that PSF, DCD, etc contain file names. If you don't have
trajectories at hand you can play with the ones that come with MDAnalysis for
testing (see below under `Examples`_)::

  >>> u = MDAnalysis.Universe(PSF, DCD)

.. SeeAlso:: :class:`MDAnalysis.core.AtomGroup.Universe` for details


"""

__version__ = "0.1.0-dev"  # NOTE: keep in sync with RELEASE in setup.py

__all__ = ['Sim', 'SimSet', 'SimAnalysis', 'SimSetAnalysis']

# Bring some often used objects into the current namespace
from Sim import Timeseries
from core.AtomGroup import Universe, asUniverse
from coordinates.core import writer as Writer

