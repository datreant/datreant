"""
AggLimbs are limbs specifically built for collections, in particular
:class:`datreant.core.collections.Bundle`. They often correspond to Treant
limbs but serve as aggregators over collections of them.

"""
import six

from . import filesystem
from . import collections
from . import _AGGLIMBS


class _AggLimbmeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)

        limbname = classdict['_name']
        _AGGLIMBS[limbname] = cls


class AggLimb(six.with_metaclass(_AggLimbmeta, object)):
    """Core functionality for limbs attached to a collection.

    """
    _name = 'agglimb'

    def __init__(self, members):
        self._members = members
