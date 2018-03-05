"""Tests for Trees and Leaves.

"""

import pytest
import os
import py

import datreant as dtr
from datreant import Veg, Leaf, Tree, Treant


class TestVeg(object):
    """Common element tests of Trees and Leaves"""
    cls = Veg
    name = 'veggie'

    @pytest.fixture
    def veg(self, tmpdir):
        with tmpdir.as_cwd():
            v = Veg(self.name)
            yield v

    def test_str(self, veg):
        assert str(veg) == os.path.join(os.getcwd(), self.name)

    def test_hash(self, veg):
        # hases are based only on abspath
        vset = {veg}
        assert veg in vset
        v2 = self.cls(self.name)
        assert v2 in vset

    def test_exists(self, veg):
        assert not veg.exists

    def test_abspath(self, veg):
        assert veg.abspath == os.path.join(os.getcwd(), self.name)

    def test_relpath(self, veg):
        assert veg.relpath == self.name

    def test_parent(self, veg):
        p = veg.parent
        assert isinstance(p, Tree)
        assert p == Tree(os.path.split(self.name)[0])

    def test_name(self, veg):
        assert veg.name == os.path.split(self.name)[1]


class TestTree(TestVeg):
    """Test generic Treant features"""
    cls = Tree
    name = 'testtreant'

    @pytest.fixture
    def tree(self, tmpdir):
        with tmpdir.as_cwd():
            t = Tree(self.name)
            yield t
    veg = tree

    def test_abspath(self, veg):
        assert veg.abspath == os.path.join(os.getcwd(), self.name) + os.sep

    def test_relpath(self, veg):
        assert veg.relpath == self.name + os.sep

    class TestTreeInit(object):
        """Test tree init.

        Test that tree works for:
            1. nonexistent directory
            2. existing directory

        Test that exception raised for:
            1. tree initialized with existing file

        """
        def test_nonexistant(self, tmpdir):
            with tmpdir.as_cwd():
                # test nonexistent directory
                t = Tree('bark')
                assert isinstance(t, Tree)
                assert not t.exists

        def test_existing(self, tmpdir):
            with tmpdir.as_cwd():
                os.makedirs('bark/lark')
                t = Tree('bark/lark/')
                assert isinstance(t, Tree)
                assert t.exists

        def test_file_ValueError(self, tmpdir):
            with tmpdir.as_cwd():
                os.mkdir('bark')
                with open(os.path.join('bark', 'mark.txt'), 'w') as f:
                    f.write('hello\nthis is a cool file\n')
                with pytest.raises(ValueError):
                    t = Tree('bark/mark.txt')

    def test_exists(self, tree):
        tree.make()

        assert os.path.exists(tree.abspath)
        assert tree.exists is True

    def test_not_exists(self, tree):
        assert not tree.exists

    @pytest.fixture
    def contains_Tree(self, tmpdir):
        # Contains various combinations of files and directories
        # that do and do not exist. Structure:
        # container
        # + dir1
        #   + file2
        # + dir3  # doesn't exist
        #   + file4  # doesn't exist
        with tmpdir.as_cwd():
            os.mkdir('container')
            os.mkdir(os.path.join('container', 'dir1'))
            with open(os.path.join('container', 'dir1', 'file2'), 'w') as f:
                f.write('some data here\n')
            tree = Tree('container')
            yield tree

    @pytest.mark.parametrize('path,exp', (
        ('dir1/', True),
        ('dir1/file2', True),
        ('dir3/', False),
        ('dir3/file4', False),
    ))
    def test_exists(self, contains_Tree, path, exp):
        assert contains_Tree[path].exists == exp

    @pytest.mark.parametrize('path,exp', (
        (os.path.join(name, 'thing1'), True),
        (os.path.join(name, 'thing2', 'thing3'), True),
        (os.path.join('other', 'thing1'), False),
    ))
    # loop over possible input types
    @pytest.mark.parametrize('inptype', (Leaf, Tree, str))
    def test_contains(self, tree, inptype, path, exp):
        thing = inptype(path)  # convert to desired type
        assert (thing in tree) == exp

    def test_contains_TypeError(self, tree):
        with pytest.raises(TypeError):
            24.0 in tree

    def test_loc(self, tree):
        l = tree.loc
        assert isinstance(l, dtr.trees._Loc)

    @pytest.fixture(params=[
        lambda x: x,  # ie just pass on the tree
        lambda x: getattr(x, 'loc')  # pass on the .loc of it
    ])
    def tree_getitem(self, tree, request):
        # for testing getitem, yields either the Tree or the '.loc' of it
        # these should be equivalent
        yield request.param(tree)

    def test_getitem_subtree(self, tree, tree_getitem):
        subt = tree_getitem['ground/control/to/major/treebeard/']
        assert isinstance(subt, Tree)
        assert not subt.exists
        assert subt.path
        assert subt in tree

    def test_getitem_leaf(self, tree, tree_getitem):
        leaf = tree_getitem['this/is/a/file']
        assert isinstance(leaf, Leaf)
        assert leaf in tree

    def test_getitem_many_leaves(self, tree_getitem):
        v = tree_getitem[['a/file', 'a/tree/']]
        assert len(v) == 2
        assert len(v.memberleaves) == 1
        assert len(v.membertrees) == 1

    def test_getitem_ValueError(self, tree_getitem):
        with pytest.raises(ValueError):
            tree_getitem['lolcats', 'a/not/file']

    def test_getitem_returntypes(self, tree):
        tree['ground/hogs/on/mars/'].make()
        tree['the/file/of your/life'].make()
        v = tree[['ground/hogs/on/mars', 'the/file/of your/life']]

        assert isinstance(v[0], Tree)
        assert isinstance(v[1], Leaf)

    def test_leaves(self, tree):
        with pytest.raises(OSError):
            tree.leaves()

        # actually make the directory now
        tree.makedirs()

        tree['.hide/me'].make()
        tree['.hide/here/'].make()

        assert len(tree.leaves()) == 0

        tree['thing1'].make()
        tree['thing2'].make()
        tree['thing3'].make()

        assert len(tree.leaves()) == 3

        tree['larry/'].make()
        tree['curly/'].make()

        assert len(tree.leaves()) == 3

    def test_trees(self, tree):
        with pytest.raises(OSError):
            tree.trees()

        # actually make the directory now
        tree.makedirs()

        assert len(tree.trees()) == 0

        tree['thing1'].make()
        tree['thing2'].make()
        tree['thing3'].make()

        assert len(tree.trees()) == 0

        tree['larry/'].make()
        tree['curly/'].make()

        assert len(tree.trees()) == 2

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

    def test_glob_OSError(self, tree):
        with pytest.raises(OSError) as error:
            tree.glob('something.*')
        assert "Tree doesn't exist in the filesystem" in str(error)

    def test_walk(self, tree):
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

    def test_walk_OSError(self, tree):
        with pytest.raises(OSError) as error:
            for v in tree.walk():  # need to use generator to trigger OSError?!
                assert v == 1
        assert "Tree doesn't exist in the filesystem" in str(error)

    def test_draw_OSError(self, tree):
        with pytest.raises(OSError) as error:
            tree.draw()
        assert "Tree doesn't exist in the filesystem" in str(error)


class TestLeaf(TestVeg):
    """Test Leaf-specific features.

    """
    cls = Leaf
    name = 'treeroot/a/fault/in/our/roots'

    @pytest.fixture
    def leaf(self, tmpdir):
        with tmpdir.as_cwd():
            l = Leaf(self.name)
            yield l
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


@pytest.fixture
def draw_tree(tmpdir):
    # set up a tree to draw with
    with tmpdir.as_cwd():
        t = Tree('here')
        t['file_zero'].make()
        t['.file_zero_hidden'].make()
        t['dir_one/'].make()
        t['dir_one/file_one'].make()
        t['dir_one/.file_one_hidden'].make()
        t['dir_one/dir_two/'].make()
        t['dir_one/dir_two/file_two'].make()
        t['dir_one/dir_two/.file_two_hidden'].make()

        yield t

DRAWREF_d0_T = """\
here/
 +-- .file_zero_hidden
 +-- file_zero
 +-- dir_one/
     +-- .file_one_hidden
     +-- file_one
     +-- dir_two/
         +-- .file_two_hidden
         +-- file_two
"""

DRAWREF_d0_F = """\
here/
 +-- file_zero
 +-- dir_one/
     +-- file_one
     +-- dir_two/
         +-- file_two
"""

DRAWREF_d1_T = """\
here/
 +-- .file_zero_hidden
 +-- file_zero
 +-- dir_one/
"""

DRAWREF_d1_F = """\
here/
 +-- file_zero
 +-- dir_one/
"""

DRAWREF_d2_T = """\
here/
 +-- .file_zero_hidden
 +-- file_zero
 +-- dir_one/
     +-- .file_one_hidden
     +-- file_one
     +-- dir_two/
"""

DRAWREF_d2_F = """\
here/
 +-- file_zero
 +-- dir_one/
     +-- file_one
     +-- dir_two/
"""


@pytest.mark.parametrize('depth,hidden,ref', (
    (None, True, DRAWREF_d0_T),
    (None, False, DRAWREF_d0_F),
    (1, True, DRAWREF_d1_T),
    (1, False, DRAWREF_d1_F),
    (2, True, DRAWREF_d2_T),
    (2, False, DRAWREF_d2_F),
    (3, True, DRAWREF_d0_T),  # 3 is max depth, so goes to full
    (3, False, DRAWREF_d0_F),
))
def test_tree_draw(draw_tree, depth, hidden, ref, capsys):
    draw_tree.draw(depth=depth, hidden=hidden)
    out, err = capsys.readouterr()
    assert out == ref
