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


def test_discover_depth(tmpdir):
    """Check that using `depth` parameter gives expected result."""
    with tmpdir.as_cwd():

        ghosts = ('something/inky',
                  'something/else/blinky',
                  'pinky',
                  'something/clyde')

        for name in ghosts:
            dtr.Treant(name)

        assert len(discover('.', depth=0)) == 0
        assert len(discover('pinky', depth=0)) == 1

        assert len(discover('.', depth=1)) == 1
        assert len(discover('.', depth=2)) == 3
        assert len(discover('.', depth=3)) == 4


def test_discover_treantdepth(tmpdir):
    """Check that using `treantdepth` parameter gives expected result."""
    with tmpdir.as_cwd():

        ghosts = ('inky',
                  'inky/blinky',
                  'pinky',
                  'inky/blinky/clyde')

        for name in ghosts:
            dtr.Treant(name)

        assert len(discover('.', treantdepth=0)) == 2
        assert len(discover('pinky', treantdepth=0)) == 1
        assert len(discover('inky', treantdepth=0)) == 1

        assert len(discover('.', treantdepth=1)) == 3
        assert len(discover('.', treantdepth=2)) == 4
        assert len(discover('inky', treantdepth=1)) == 2
        assert len(discover('inky/blinky', treantdepth=1)) == 2

        assert len(discover('inky/blinky', treantdepth=1)) == 2


def test_discover_depth_treantdepth(tmpdir):
    """Check that using `treantdepth` and `depth` parameters together gives
        expected result.
    """
    with tmpdir.as_cwd():

        ghosts = ('inky',
                  'inky/blinky',
                  'pinky',
                  'inky/blinky/nothing/clyde')

        for name in ghosts:
            dtr.Treant(name)

        assert len(discover('.', treantdepth=0, depth=0)) == 0
        assert len(discover('.', treantdepth=0, depth=1)) == 2
        assert len(discover('pinky', treantdepth=0, depth=0)) == 1
        assert len(discover('inky', treantdepth=0, depth=2)) == 1

        assert len(discover('.', treantdepth=1, depth=1)) == 2
        assert len(discover('inky', treantdepth=1, depth=1)) == 2
        assert len(discover('inky', treantdepth=1, depth=0)) == 1

        assert len(discover('inky', treantdepth=2)) == 3
        assert len(discover('inky', treantdepth=2, depth=2)) == 2
        assert len(discover('inky', treantdepth=2, depth=3)) == 3
