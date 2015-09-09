# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# datreant

"""
datreant --- a persistence engine for heterogeneous data sets
=============================================================

.. SeeAlso:: :class:`datreant.treants.Treant`
             :class:`datreant.treants.Group`

"""
# Bring some often used objects into the current namespace
from datreant.treants import Treant, Group, register
from datreant.coordinator import Coordinator
from datreant.collections import Bundle
from datreant.manipulators import *
import datreant.persistence

__all__ = ['Treant', 'Group', 'Coordinator', 'Bundle']
__version__ = "0.6.0-dev"  # NOTE: keep in sync with RELEASE in setup.py

_treants = dict()
register(Treant, Group)
