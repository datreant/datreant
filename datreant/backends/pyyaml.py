"""
Interface classes for YAML state files.

"""

import os
import sys
import fcntl
import pickle
import logging
import warnings
from functools import wraps

import yaml

import datreant
from .serial import TreantFileSerial, GroupFileSerial


class MixinYAML(object):

    def _deserialize(self, handle):
        return yaml.load(handle)

    def _serialize(self, record, handle):
        yaml.dump(record, handle)


class TreantFileYAML(MixinYAML, TreantFileSerial):
    pass


class GroupFileYAML(MixinYAML, GroupFileSerial):
    pass
