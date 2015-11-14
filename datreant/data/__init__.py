# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# datreant

"""
datreant.data --- basic data storage backends for datreant.limbs.Data
=====================================================================

"""
from datreant.data.core import DataFile
from datreant.data import pydata, npdata, pddata

__all__ = ['pydata', 'npdata', 'pddata']
