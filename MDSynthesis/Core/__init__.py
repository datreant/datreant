# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDSynthesis.Core

"""
:mod:`MDSynthesis.Core` --- Lower-level, non-user classes
==========================================================

This submodule contains all the elements that support the Containers, Operators,
and Databases of MDSynthesis. These classes are not intended for direct use
by users, but some may be used as interfaces in user-level objects.

"""
# Bring everything into the Core namespace
import Aggregators 
import Files
import Workers
