"""Tests for Trees and Leaves.

"""

import pytest
import os
import py

from datreant.core import Veg, Leaf, Tree, Treant


class TestVeg:
    """Common element tests of Trees and Leaves"""

    @pytest.fixture
    def veg(self, tmpdir):
        with tmpdir.as_cwd():
            v = Veg('veggie')
        return v


class TestTree(TestVeg):
    """Test generic Treant features"""
    treeroot = 'testtreant'

    @pytest.fixture
    def tree(self, tmpdir):
        with tmpdir.as_cwd():
            t = Tree(self.treeroot)
        return t

    veg = tree

    def test_init(self, tmpdir):
        """Test tree init.

        Test that tree works for:
            1. nonexistent directory
            2. existing directory

        Test that exception raised for:
            1. tree initialized with existing file

        """
        with tmpdir.as_cwd():

            # test nonexistent directory
            t = Tree('bark')
            assert not t.exists

            # test existent directory
            t2 = t['lark/'].make()
            assert t2.exists

            t3 = Tree('bark/lark')
            assert t3.exists

            # test that init with file raises ValueError

            # this should create a nonexistent Tree
            t4 = Tree('bark/mark.txt')
            assert not t4.exists

            # this makes a file
            t['mark.txt'].make()

            with pytest.raises(ValueError):
                t5 = Tree('bark/mark.txt')

    def test_exists(self, tree):
        tree.make()

        assert os.path.exists(tree.abspath)
        assert tree.exists is True

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

        v = tree[['a/file', 'a/tree/']]

        assert len(v) == 2
        assert len(v.memberleaves) == 1
        assert len(v.membertrees) == 1

        with pytest.raises(ValueError):
            tree['lolcats', 'a/not/file']

        tree['ground/hogs/on/mars/'].make()
        tree['the/file/of your/life'].make()

        v = tree[['ground/hogs/on/mars', 'the/file/of your/life']]

        assert isinstance(v[0], Tree)
        assert isinstance(v[1], Leaf)

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

    def test_walk(self, tree, tmpdir):

        with tmpdir.as_cwd():
            # files
            tree['scipy'].make()
            tree['2016'].make()
            tree['sprint'].make()

            # directories with files
            t1 = tree['a_dir/has_no_name']
            t2 = tree['another_dir/bites_the_dust']
            t1.make()
            t2.make()

            roots_scandir = []
            dirs_scandir = []
            files_scandir = []
            all_roots = []
            all_trees = []
            all_leaves = []

            for root, dirs, files in os.walk(tree.abspath):
                # os.walk normally doesn't add slash to path, in order to
                # replicate the path given by tree.walk() we use os.path.join
                if root != tree.abspath:
                    root = os.path.join(root, '')
                roots_scandir.append(root)
                for directory in dirs:
                    # this is the abspath of the directory,
                    # same reason as above for use of os.path.join
                    dirs_scandir.append(os.path.join(root, directory, ''))
                for f in files:
                    files_scandir.append(f)

            for root, trees, leaves in tree.walk():
                all_roots.append(root.abspath)
                for tree in trees:
                    # this is the abspath of the directory
                    all_trees.append(tree.abspath)
                for leaf in leaves:
                    all_leaves.append(leaf.name)

            assert roots_scandir == all_roots
            assert dirs_scandir == all_trees
            assert files_scandir == all_leaves


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

    def test_init(self, tmpdir):
        """
        Test that leaf works for:
            1. nonexistent file
            2. existing file

        Test that exception raised for:
            1. leaf initialized with existing directory

        """
        with tmpdir.as_cwd():

            # test nonexistent file
            t = Leaf('bark')
            assert not t.exists

            # test existent file
            t.make()
            assert t.exists

            t2 = Leaf('bark')
            assert t2.exists

            # test that init with directory raises ValueError

            # this should create a nonexistent Tree
            t3 = Tree('mark/').make()
            assert t3.exists

            with pytest.raises(ValueError):
                t4 = Leaf('mark')

    def test_exists(self, leaf):
        leaf.make()

        assert os.path.exists(leaf.abspath)
        assert leaf.exists is True

    def test_makedirs(self, leaf):
        pass
