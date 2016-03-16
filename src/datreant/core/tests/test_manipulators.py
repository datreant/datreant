"""Tests for manipulators.

"""

import pytest

import datreant.core as dtr
from datreant.core.manipulators import discover


def test_discover(tmpdir):
    with tmpdir.as_cwd():

        ghosts = ('inky', 'blinky', 'pinky', 'clyde')

        for name in ghosts:
            dtr.Treant(
                    'a/very/deep/directory/structure/that/just/keeps/going/' +
                    name)

        b = discover('.')

        assert len(b) == 4

        for name in ghosts:
            assert name in b.names
