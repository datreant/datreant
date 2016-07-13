# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# datreant.core

"""
datreant.core --- persistent, pythonic trees for heterogeneous data
===================================================================

.. SeeAlso:: :class:`datreant.core.treants.Treant`
             :class:`datreant.core.treants.Group`

"""
# global registries of classes
_TREANTS = dict()

_TREELIMBS = dict()
_LIMBS = dict()
_AGGTREELIMBS = dict()
_AGGLIMBS = dict()

# Bring some often used objects into the current namespace
from .manipulators import discover
from .treants import Treant, Group
from .trees import Veg, Leaf, Tree
from .collections import View, Bundle
from . import attach

__all__ = ['Treant', 'Group', 'Tree', 'Leaf', 'Bundle']
__version__ = "0.7.0"  # NOTE: keep in sync with RELEASE in setup.py
