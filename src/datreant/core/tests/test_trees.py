"""Tests for Trees and Leaves.

"""

import pytest
import os
import py
from . import test_bundle

from datreant.core import Veg, Leaf, Tree


class TestVeg:
    """Common element tests of Trees and Leaves"""

    def test_exists(self, veg):
        pass


class TestTree(TestVeg):
    """Test generic Treant features"""
    treeroot = 'testtreant'

    @pytest.fixture
    def tree(self, tmpdir):
        with tmpdir.as_cwd():
            t = Tree(self.treeroot)
        return t

    veg = tree

    def test_getitem(self, tree):
        """Test that using getitem syntax returns Trees and Leaves as it
        should.

        """
        subt = tree['ground/control/to/major/treebeard/']
        assert isinstance(subt, Tree)
        assert not subt.exists
        assert subt in tree

        leaf = tree['this/is/a/file']
        assert isinstance(leaf, Leaf)
        assert leaf in tree

    def test_leaves(self, tree):
        with pytest.raises(OSError):
            tree.leaves

        # actually make the directory now
        tree.makedirs()

        tree['.hide/me'].make()
        tree['.hide/here/'].make()

        assert len(tree.leaves) == 0

        tree['thing1'].make()
        tree['thing2'].make()
        tree['thing3'].make()

        assert len(tree.leaves) == 3

        tree['larry/'].make()
        tree['curly/'].make()

        assert len(tree.leaves) == 3

    def test_trees(self, tree):
        with pytest.raises(OSError):
            tree.trees

        # actually make the directory now
        tree.makedirs()

        assert len(tree.trees) == 0

        tree['thing1'].make()
        tree['thing2'].make()
        tree['thing3'].make()

        assert len(tree.trees) == 0

        tree['larry/'].make()
        tree['curly/'].make()

        assert len(tree.trees) == 2

    def test_hidden(self, tree):
        with pytest.raises(OSError):
            tree.hidden

        # actually make the directory now
        tree.makedirs()

        assert len(tree.hidden) == 0

        tree['.hide/me'].make()
        tree['.hide/here/'].make()
        tree['.hide.nothing/here/'].make()

        assert len(tree.hidden) == 2

        tree['thing1'].make()
        tree['thing2'].make()
        tree['thing3'].make()

        assert len(tree.hidden) == 2

        tree['larry/'].make()
        tree['curly/'].make()

        assert len(tree.hidden) == 2

    # def test_draw(self, tree):
    #     from io import StringIO
    #     import sys

    #     with pytest.raises(OSError):
    #         tree.draw()

    #     class Capturing(list):
    #         def __enter__(self):
    #             self._stdout = sys.stdout
    #             sys.stdout = self._stringio = StringIO()
    #             return self
    #         def __exit__(self, *args):
    #             self.extend(self._stringio.getvalue().splitlines())
    #             sys.stdout = self._stdout

    #     tree.makedirs()

    #     with Capturing() as output:
    #         tree.draw()

    #     assert tree.relpath in output[0]

    def test_makedirs(self, tree):
        pass


class TestLeaf(TestVeg):
    """Test Leaf-specific features.

    """
    leafname = 'treeroot/a/fault/in/our/roots'

    @pytest.fixture
    def leaf(self, tmpdir):
        with tmpdir.as_cwd():
            l = Leaf(self.leafname)
        return l

    veg = leaf

    def test_makedirs(self, leaf):
        pass
