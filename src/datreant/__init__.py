# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# datreant

"""
datreant --- persistent, pythonic trees for heterogeneous data
==============================================================

.. SeeAlso:: :class:`datreant.treants.Treant`

"""

# Bring some often used objects into the current namespace
from .manipulators import discover
from .treants import Treant
from .trees import Veg, Leaf, Tree
from .collections import View, Bundle

__all__ = ['Treant', 'Tree', 'Leaf', 'Bundle', 'discover', 'View']
__version__ = "1.0.2"  # NOTE: keep in sync with RELEASE in setup.py
