"""
Interface classes for JSON state files.

"""

import os
import sys
import fcntl
import pickle
import logging
import warnings
from functools import wraps

import json

import datreant
from .serial import TreantFileSerial


class TreantFileJSON(TreantFileSerial):

    def _deserialize(self, handle):
        return json.load(handle)

    def _serialize(self, record, handle):
        json.dump(record, handle)
