"""Tests for Trees and Leaves.

"""

import pytest
import os
import py
from . import test_bundle

from datreant.core import Tree, Leaf


class Brush:
    """Common element tests of Trees and Leaves"""

    def test_exists(self, brush):
        brush.make()

        assert os.path.exists(brush.abspath)
        assert brush.exists is True


class TestTree(Brush):
    """Test generic Treant features"""
    treeroot = 'testtreant'

    @pytest.fixture
    def tree(self, tmpdir):
        with tmpdir.as_cwd():
            t = Tree(self.treeroot)
        return t

    brush = tree

    def test_init(self, tmpdir):
        """
        Test that tree works for:
            1. nonexistent directory
            2. existing directory

        Test that exception raised for:
            1. tree initialized with existing file

        """
        pass

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

    def test_equal(self, tree):
        t1 = tree['a dir/']
        t2 = tree['another dir/']

        assert t1 != t2
        assert t1['.'] == t1

    def test_compare(self, tree):
        assert tree['bark/'] <= tree['dark/']

    def test_makedirs(self, tree):
        t1 = tree['a/ton of/stupid/bricks/'].makedirs()

        assert t1.exists

    def test_glob(self, tree):
        tm = tree['moe'].make()
        tl = tree['larry'].make()
        tc = tree['curly'].make()

        assert tl in tree.glob('*r*y')
        assert tc in tree.glob('*r*y')


class TestLeaf(Brush):
    """Test Leaf-specific features.

    """
    leafname = 'treeroot/a/fault/in/our/roots'

    @pytest.fixture
    def leaf(self, tmpdir):
        with tmpdir.as_cwd():
            l = Leaf(self.leafname)
        return l

    brush = leaf

    def test_init(self, tmpdir):
        """
        Test that leaf works for:
            1. nonexistent file
            2. existing file

        Test that exception raised for:
            1. leaf initialized with existing directory

        """
        pass

    def test_makedirs(self, leaf):
        pass
