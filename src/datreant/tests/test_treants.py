"""Interface tests for Treants.

"""

import datreant as dtr
from datreant import Treant
import pytest
import mock
import os
import py

from . import test_collections
from .test_trees import TestTree


class TestTreant(TestTree):
    """Test generic Treant features"""
    treantname = 'testtreant'

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.treants.Treant(TestTreant.treantname)
        return c

    @pytest.fixture
    def basic_treant(self, tmpdir):
        # treant with tags and cats, in tmpdir
        with tmpdir.as_cwd():
            t1 = Treant('Rincewind')
            t1.tags.add('magical')
            t1.categories.add({'colour': 'octarine'})
            yield t1

    def test_init(self, treant, tmpdir):
        """Test basic Treant init"""
        assert treant.name == self.treantname
        assert treant.abspath == (os.path.join(tmpdir.strpath,
                                               self.treantname) + os.sep)

    def test_init_from_Tree(self, tmpdir):
        with tmpdir.as_cwd():
            tree = dtr.Tree('this')
            t = Treant(tree)

            assert t.path == tree.path

    @pytest.mark.parametrize("tags", (None, [], ['small', 'spiky']))
    @pytest.mark.parametrize("categories", (None, {}, {'colour': 'red'}))
    def test_init_generate(self, tags, categories, tmpdir):
        # test combinations of tags and categories
        # when generating from scratch
        with tmpdir.as_cwd():
            t = dtr.Treant('babs', tags=tags, categories=categories)

            if tags is not None:
                for tag in tags:
                    assert tag in t.tags
            if categories is not None:
                for cat, val in categories.items():
                    assert cat in t.categories
                    assert t.categories[cat] == val

    @pytest.mark.parametrize("tags", (None, [], ['small', 'spiky']))
    @pytest.mark.parametrize("categories", (None, {}, {'colour': 'red'}))
    def test_init_regenerate_via_name(self, tags, categories, tmpdir):
        # test regenerating a Treant from its directory
        with tmpdir.as_cwd():
            t = dtr.Treant('this')

            t2 = dtr.Treant('this', tags=tags, categories=categories)
            if tags is not None:
                for tag in tags:
                    assert tag in t2.tags
            if categories is not None:
                for cat, val in categories.items():
                    assert cat in t2.categories
                    assert t2.categories[cat] == val

    def test_gen_OSError(self, tmpdir):
        with tmpdir.as_cwd():
            with mock.patch('os.makedirs') as mp:
                mp.sideeffect = OSError(os.errno.ENOSPC, 'Mock - disk full')
                with pytest.raises(OSError) as error:
                    t = Treant('new')
                    t.tags.add('worthless')
                    assert error.errno == os.errno.ENOSPC

    def test_gen_OSError13(self, tmpdir):
        with tmpdir.as_cwd():
            with mock.patch('os.makedirs') as mp:
                mp.sideeffect = OSError(os.errno.EACCES, 'Mock - disk full')
                with pytest.raises(OSError) as error:
                    t = Treant('new')
                    t.tags.add('worthless')
                    assert error.errno == os.errno.EACCES
                    assert ("Permission denied; cannot create 'new'"
                            in str(error))

    def test_gen_methods(self, tmpdir):
        """Test the variety of ways we can generate a new Treant

        1. ``Treant('treant')``, where 'treant' is not an existing file or
           directory path

        2. ``Treant('treant')``, where 'treant' is an existing directory
           without Treant state files inside

        3. ``Treant('/somedir/treant')``, where 'treant' is not an existing
           file or directory in 'somedir'

        4. ``Treant('/somedir/treant')``, where 'treant' is an existing
           directory in 'somedir' without any Treant state files inside

        """
        with tmpdir.as_cwd():
            # 1
            t1 = Treant('newone')
            assert os.path.exists(t1._treantdir)

            # 2
            os.mkdir('another')
            t2 = Treant('another')
            assert os.path.exists(t2._treantdir)

            # 3
            t3 = Treant('yet/another')
            assert os.path.exists(t3._treantdir)

            # 4
            os.mkdir('yet/more')
            t4 = Treant('yet/more')
            assert os.path.exists(t4._treantdir)

    @pytest.mark.parametrize("tags", (None, [], ['small', 'spiky']))
    @pytest.mark.parametrize("categories", (None, {}, {'colour': 'red'}))
    def test_regen(self, tags, categories, tmpdir):
        """Test regenerating Treant.

        - create Treant
        - modify Treant a little
        - create same Treant (should regenerate)
        - check that modifications were saved
        """
        with tmpdir.as_cwd():
            C1 = Treant('regen', tags=tags, categories=categories)

            C2 = Treant('regen')  # should be regen of C1

            if tags is not None:
                for tag in tags:
                    assert tag in C2.tags
            if categories is not None:
                for cat, val in categories.items():
                    assert cat in C2.categories
                    assert C2.categories[cat] == val

            assert C1 == C2

            # they point to the same file, but they are not the same object
            assert C1 is not C2

    def test_cmp(self, tmpdir):
        """Test the comparison of Treants when sorting"""
        with tmpdir.as_cwd():
            c1 = Treant('a')
            c2 = Treant('b')
            c3 = Treant('c')

        assert sorted([c3, c2, c1]) == [c1, c2, c3]
        assert c1 <= c2 < c3
        assert c3 >= c2 > c1

    class TestTags:
        """Test treant tags"""

        @pytest.mark.parametrize("tags", (("1", "2", "3"),   # tuple
                                          ["1", "2", "3"],   # list
                                          {"1", "2", "3"}))  # set
        def test_adding_array_like_tags(self, treant, tags):
            treant.tags.add(tags)
            assert sorted(list(tags)) == sorted(list(treant.tags))

        def test_add_tags(self, treant):
            treant.tags.add('marklar')
            assert 'marklar' in treant.tags

            treant.tags.add('lark', 'bark')
            assert 'marklar' in treant.tags
            assert 'lark' in treant.tags
            assert 'bark' in treant.tags

        def test_remove_tags(self, treant):
            treant.tags.add('marklar')
            assert 'marklar' in treant.tags
            treant.tags.remove('marklar')
            assert 'marklar' not in treant.tags

            treant.tags.add('marklar')
            treant.tags.add('lark', 'bark')
            treant.tags.add(['fark', 'bark'])
            assert 'marklar' in treant.tags
            assert 'lark' in treant.tags
            assert 'bark' in treant.tags
            assert 'fark' in treant.tags
            assert len(treant.tags) == 4

            treant.tags.remove('fark')
            assert 'fark' not in treant.tags
            assert len(treant.tags) == 3
            treant.tags.remove('fark')
            assert len(treant.tags) == 3

            treant.tags.clear()
            assert len(treant.tags) == 0

        def test_tags_set_behavior(self, tmpdir):
            with tmpdir.as_cwd():
                # 1
                t1 = Treant('maple')
                t1.tags.add(['sprout', 'deciduous'])

                # 2
                t2 = Treant('sequoia')
                t2.tags.add(['sprout', 'evergreen'])

                tags_union = t1.tags | t2.tags
                for t in ['sprout', 'deciduous', 'evergreen']:
                    assert t in tags_union

                tags_intersect = t1.tags & t2.tags
                assert 'sprout' in tags_intersect
                for t in ['deciduous', 'evergreen']:
                    assert t not in tags_intersect

                tags_diff = t1.tags - t2.tags
                assert 'deciduous' in tags_diff
                for t in ['sprout', 'evergreen']:
                    assert t not in tags_diff

                tags_symm_diff = t1.tags ^ t2.tags
                for t in ['deciduous', 'evergreen']:
                    assert t in tags_symm_diff
                assert 'sprout' not in tags_symm_diff

                # 3
                t3 = Treant('oak')
                t3.tags.add(['deciduous'])

                # Test set membership
                assert t1.tags <= t1.tags
                assert not t1.tags < t1.tags
                assert t1.tags == t1.tags
                assert not t1.tags < t3.tags
                assert t1.tags > t3.tags

                # test TypeErrors in Tags
                # type_error_msg = "Operands must be AggTags, Tags, or a set."
                with pytest.raises(TypeError) as e:
                    ('tree') == t1.tags
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    t1.tags < ('tree')
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    ('tree') - t1.tags
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    t1.tags - ('tree')
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    ('tree') | t1.tags
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    t1.tags | ('tree')
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    ('tree') & t1.tags
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    t1.tags & ('tree')
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    ('tree') ^ t1.tags
                # assert e.value.message == type_error_msg
                with pytest.raises(TypeError) as e:
                    t1.tags ^ ('tree')
                # assert e.value.message == type_error_msg

        def test_tags_setting(self, tmpdir):
            """Test that we can set tags with lists or sets, or with Tags
            objects.

            """
            with tmpdir.as_cwd():

                # set with a list
                t1 = Treant('maple')
                t1.tags = ['sprout', 'deciduous']

                assert t1.tags == {'sprout', 'deciduous'}

                # set with a set
                t2 = Treant('elm')
                t2.tags = {'sprout', 'deciduous'}

                assert t2.tags == {'sprout', 'deciduous'}

                # set with a Tags object
                t3 = Treant('sequoia')
                t3.tags = t2.tags

                assert t3.tags == {'sprout', 'deciduous'}

        def test_tags_fuzzy(self, tmpdir, treant):
            """Test that fuzzy matching for tags works as expected.
            """
            treant.tags.add('bark', 'leafy', 'green', 'Leafy')

            for tag in ('leafy', 'Leafy'):
                assert tag in treant.tags.fuzzy('leafy')

        def test_tags_getitem(self, treant):
            """Test the queryability of tags via its __getitem__ method."""
            treant.tags.add('marklar', 'lark', 'bark')

            t = treant

            # single presence
            assert t.tags['lark']
            assert not t.tags['mark']

            # single not present
            assert t.tags[{'mark'}]

            # anding
            assert t.tags[['marklar', 'bark']]

            # oring
            assert t.tags['marklar', 'bark']
            assert t.tags['mark', 'bark']
            assert not t.tags['mark', 'dark']

            # not anding
            assert t.tags[{'dark', 'marklar'}]

            # complex logic
            assert t.tags[[('marklar', 'bark'), {'dark'}]]

        @pytest.mark.parametrize('tag', (1, 1.2))
        def test_tags_only_strings(self, treant, tag):
            with pytest.raises(ValueError):
                treant.tags.add(tag)

        def test_tags_printing(self, treant):
            treant.tags.add('marklar')
            repr = str(treant.tags)
            assert repr == "Tags\n====\n'marklar'\n"

    class TestCategories:
        """Test treant categories"""

        def test_add_categories(self, treant):
            treant.categories.add(marklar=42)
            assert 'marklar' in treant.categories

            treant.categories.add({'bark': 'snark'}, lark=27)
            assert 'bark' in treant.categories
            assert 'snark' not in treant.categories
            assert 'bark' in treant.categories

            assert treant.categories['bark'] == 'snark'
            assert treant.categories['lark'] == 27

            treant.categories['lark'] = 42
            assert treant.categories['lark'] == 42

        def test_remove_categories(self, treant):
            treant.categories.add(marklar=42)
            assert 'marklar' in treant.categories

            treant.categories.remove('marklar')
            assert 'marklar' not in treant.categories

            treant.categories.add({'bark': 'snark'}, lark=27)
            del treant.categories['bark']
            assert 'bark' not in treant.categories

            # should just work, even if key isn't present
            treant.categories.remove('smark')

            treant.categories['lark'] = 42
            treant.categories['fark'] = 32.3

            treant.categories.clear()
            assert len(treant.categories) == 0

        def test_add_wrong(self, treant):
            with pytest.raises(TypeError):
                treant.categories.add('temperature', 300)

            with pytest.raises(TypeError):
                treant.categories.add(['mark', 'matt'])

        @pytest.mark.parametrize('key, val', [[2, 'twenty'],
                                              [['blarg'], "nothin'"],
                                              [None, "literally nothin'"],
                                              [True, 'tautologically']])
        def test_add_wrong_keys(self, treant, key, val):
            # always test both addition methods
            with pytest.raises(TypeError):
                treant.categories[key] = val
            with pytest.raises(TypeError):
                treant.categories.add(key, val)

        @pytest.mark.parametrize('key, val', [['bark', ['shaggy']],
                                              ['snark', {'yes'}]])
        def test_add_wrong_values(self, treant, key, val):
            # always test both addition methods
            with pytest.raises(TypeError):
                treant.categories.add(key, val)
            with pytest.raises(TypeError):
                treant.categories[key] = val

        def test_None_deletes(self, treant):
            """Setting a category to ``None`` should delete the value.
            """
            treant.categories['size'] = 'large'
            treant.categories['bark'] = 'smooth'
            treant.categories['bark'] = None

            assert 'bark' not in treant.categories
            assert 'size' in treant.categories

            treant.categories['bark'] = 'smooth'
            treant.categories.add(bark=None)

            assert 'bark' not in treant.categories
            assert 'size' in treant.categories

        def test_setting_to_None_VE(self, treant):
            with pytest.raises(ValueError) as err:
                treant.categories['colour'] = None
            assert "Cannot set to 'None'" in str(err)

        def test_KeyError(self, treant):
            with pytest.raises(KeyError):
                treant.categories['hello?']

        def test_get_categories(self, treant):
            treant.categories['bark'] = 'dark'
            treant.categories['leaves'] = 'many'
            treant.categories['roots'] = 'shallow'

            # get a single category
            assert treant.categories['leaves'] == 'many'

            # get multiple categories with list
            assert treant.categories[['leaves', 'bark']] == ['many', 'dark']

            # get multiple categories with set
            assert treant.categories[{'leaves', 'bark'}] == {'leaves': 'many',
                                                             'bark': 'dark'}

        def test_set_categories(self, treant, tmpdir):
            a_dict = {'leaves': 'many', 'bark': 'dark'}

            treant.categories = a_dict
            assert treant.categories == a_dict

            a_dict.update({'roots': 'shallow'})

            assert treant.categories != a_dict
            treant.categories = a_dict
            assert a_dict == treant.categories

            # test setting from other Treant's categories
            with tmpdir.as_cwd():
                s = dtr.Treant('sprout')

            s.categories = {'shoes': False, 'shirt': False}

            assert treant.categories != s.categories

            treant.categories = s.categories
            assert treant.categories == s.categories

        def test_from_treant(self, treant, tmpdir):
            with tmpdir.as_cwd():
                dtr.Treant('sprout', categories=treant.categories)


class TestReadOnly:
    """Test Treant functionality when read-only"""

    @pytest.fixture
    def treant(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = dtr.treants.Treant('testtreant')
            c.tags.add('72')
            py.path.local(c.abspath).chmod(0o0550, rec=True)

        def fin():
            py.path.local(c.abspath).chmod(0o0770, rec=True)

        request.addfinalizer(fin)

        return c

    def test_treant_read_only(self, treant):
        """Test that a read-only Treant can be accessed, but not written to.
        """
        c = dtr.treants.Treant(treant.abspath)

        assert '72' in c.tags

        with pytest.raises(OSError):
            c.tags.add('yet another')
