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
# Bring some often used objects into the current namespace
from .treants import Treant, Group, register
from .collections import Bundle
from . import attach

__all__ = ['Treant', 'Group', 'Bundle']
__version__ = "0.6.0-dev"  # NOTE: keep in sync with RELEASE in setup.py

_TREANTS = dict()
register(Treant, Group)
