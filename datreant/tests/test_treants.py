"""Interface tests for Treants.

"""

import datreant as dtr
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py
import test_data
import test_bundle


class TestTreant:
    """Test generic Treant features"""
    treantname = 'testtreant'
    treanttype = 'Treant'
    treantclass = dtr.treants.Treant

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.treants.Treant(TestTreant.treantname)
        return c

    def test_init(self, treant, tmpdir):
        """Test basic Treant init"""
        assert treant.name == self.treantname
        assert treant.location == tmpdir.strpath
        assert treant.treanttype == self.treanttype
        assert treant.basedir == os.path.join(tmpdir.strpath, self.treantname)

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

        5. ``Treant('somedir/treant', new=True)``, where 'treant' is an
           existing directory in 'somedir' with an existing Treant statefile

        6. ``Treant('/somedir/treant')``, where 'treant' is an existing
           directory in 'somedir' with other types of Treant files inside (such
           as Group)

        """
        with tmpdir.as_cwd():
            # 1
            t1 = self.treantclass('newone')
            assert os.path.exists(t1.filepath)

            # 2
            os.mkdir('another')
            t2 = self.treantclass('another')
            assert os.path.exists(t2.filepath)

            # 3
            t3 = self.treantclass('yet/another')
            assert os.path.exists(t3.filepath)

            # 4
            os.mkdir('yet/more')
            t4 = self.treantclass('yet/more')
            assert os.path.exists(t4.filepath)

            # 5
            t5 = self.treantclass('yet/more', new=True)
            assert os.path.exists(t5.filepath)
            assert t5.basedir == t4.basedir
            assert t5.filepath != t4.filepath

            # 6
            compare = t1.uuid
            os.rename(t1.filepath,
                      t1.filepath.replace(t1.treanttype, 'Another'))
            t6 = self.treantclass('newone')
            assert t6.uuid != compare

    def test_regen(self, tmpdir):
        """Test regenerating Treant.

        - create Treant
        - modify Treant a little
        - create same Treant (should regenerate)
        - check that modifications were saved
        """
        with tmpdir.as_cwd():
            C1 = self.treantclass('regen')
            C1.tags.add('fantastic')
            C2 = self.treantclass('regen')  # should be regen of C1
            assert 'fantastic' in C2.tags

            # they point to the same file, but they are not the same object
            assert C1 is not C2

    def test_regen_methods(self, tmpdir):
        """Test the variety of ways Treants can be regenerated.

        1. ``Treant('treant')``, where there exists *only one* ``Treant`` state
           file inside 'treant'

        2. ``Treant('treant/Treant.<uuid>.<ext>')``, where there need not be
           only a single ``Treant`` state file in 'treant'

        """
        with tmpdir.as_cwd():
            t1 = self.treantclass('newone')
            t2 = self.treantclass('newone')
            assert t1.uuid == t2.uuid

            t3 = self.treantclass('newone', new=True)
            assert t3.uuid != t2.uuid

            t4 = self.treantclass(t3.filepath)
            assert t4.uuid == t3.uuid

    def test_noregen(self, tmpdir):
        """Test a variety of ways that generation of a new Treant should fail.

        1. `Treant('somedir/treant')` should raise `MultipleTreantsError` if
           more than one state file is in the given directory

        """
        with tmpdir.as_cwd():
            # 1
            t1 = self.treantclass('newone')
            t2 = self.treantclass('newone', new=True)
            assert t1.uuid != t2.uuid

            with pytest.raises(dtr.treants.MultipleTreantsError):
                t3 = self.treantclass('newone')

    def test_cmp(self, tmpdir):
        """Test the comparison of Treants when sorting"""
        with tmpdir.as_cwd():
            c1 = self.treantclass('a')
            c2 = self.treantclass('b')
            c3 = self.treantclass('c')

        assert sorted([c3, c2, c1]) == [c1, c2, c3]
        assert c1 <= c2 < c3
        assert c3 >= c2 > c1

    class TestTags:
        """Test treant tags"""

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

            treant.tags.remove(all=True)
            assert len(treant.tags) == 0

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

            treant.categories.remove(all=True)
            assert len(treant.categories) == 0

        def test_add_wrong(self, treant):
            with pytest.raises(TypeError):
                treant.categories.add('temperature', 300)

            with pytest.raises(TypeError):
                treant.categories.add(['mark', 'matt'])

        def test_KeyError(self, treant):
            with pytest.raises(KeyError):
                treant.categories['hello?']

    class TestData:
        """Test data storage and retrieval"""

        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_add_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(treant.basedir,
                                                   self.handle,
                                                   self.datafile))

            def test_remove_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(treant.basedir,
                                                   self.handle,
                                                   self.datafile))

                treant.data.remove('testdata')
                assert not os.path.exists(os.path.join(treant.basedir,
                                                       self.handle,
                                                       self.datafile))

            def test_retrieve_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                np.testing.assert_equal(treant.data.retrieve(self.handle),
                                        datastruct)
                np.testing.assert_equal(treant.data[self.handle],
                                        datastruct)

        class PandasMixin(DataMixin):
            """Mixin class for pandas tests"""
            datafile = dtr.data.pddata.pddatafile

            def test_retrieve_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                np.testing.assert_equal(
                        treant.data.retrieve(self.handle).values,
                        datastruct.values)
                np.testing.assert_equal(
                        treant.data[self.handle].values,
                        datastruct.values)

        class Test_Series(test_data.Series, PandasMixin):
            pass

        class Test_DataFrame(test_data.DataFrame, PandasMixin):
            pass

        class Test_Blank_DataFrame(test_data.Blank_DataFrame, PandasMixin):
            pass

        class Test_Wide_Blank_DataFrame(test_data.Wide_Blank_DataFrame,
                                        PandasMixin):
            pass

        class Test_Thin_Blank_DataFrame(test_data.Thin_Blank_DataFrame,
                                        PandasMixin):
            pass

        class Test_Panel(test_data.Panel, PandasMixin):
            pass

        class Test_Panel4D(test_data.Panel4D, PandasMixin):
            pass

        class NumpyMixin(DataMixin):
            """Test numpy datastructure storage and retrieval"""
            datafile = dtr.data.npdata.npdatafile

        class Test_NumpyScalar(test_data.NumpyScalar, NumpyMixin):
            pass

        class Test_Numpy1D(test_data.Numpy1D, NumpyMixin):
            pass

        class Test_Numpy2D(test_data.Numpy2D, NumpyMixin):
            pass

        class Test_Wide_Numpy2D(test_data.Wide_Numpy2D, NumpyMixin):
            pass

        class Test_Thin_Numpy2D(test_data.Thin_Numpy2D, NumpyMixin):
            pass

        class Test_Numpy3D(test_data.Numpy3D, NumpyMixin):
            pass

        class Test_Numpy4D(test_data.Numpy4D, NumpyMixin):
            pass

        class PythonMixin(DataMixin):
            """Test pandas datastructure storage and retrieval"""
            datafile = dtr.data.pydata.pydatafile

            def test_overwrite_data(self, treant, datastruct):
                treant.data[self.handle] = datastruct

                # overwrite the data with a scalar
                treant.data[self.handle] = 23
                assert treant.data[self.handle] == 23

        class Test_List(test_data.List, PythonMixin):
            pass

        class Test_Dict(test_data.Dict, PythonMixin):
            pass

        class Test_Tuple(test_data.Tuple, PythonMixin):
            pass

        class Test_Set(test_data.Set, PythonMixin):
            pass

        class Test_Dict_Mix(test_data.Dict_Mix, PythonMixin):
            pass


class TestGroup(TestTreant):
    """Test Group-specific features.

    """
    treantname = 'testgroup'
    treanttype = 'Group'
    treantclass = dtr.Group

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            g = dtr.Group(TestGroup.treantname)
        return g

    def test_repr(self, treant):
        pass

    class TestMembers(test_bundle.CollectionTests):
        """Test member functionality"""

        @pytest.fixture
        def collection(self, treant):
            return treant.members


class TestReadOnly:
    """Test Treant functionality when read-only"""

    @pytest.fixture
    def treant(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = dtr.treants.Treant('testtreant')
            c.tags.add('72')
            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    @pytest.fixture
    def group(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = dtr.Group('testgroup')
            c.members.add(dtr.Treant('lark'), dtr.Group('bark'))
            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    def test_treant_read_only(self, treant):
        """Test that a read-only Treant can be accessed, but not written to.
        """
        c = dtr.treants.Treant(treant.filepath)

        assert '72' in c.tags

        with pytest.raises(OSError):
            c.tags.add('yet another')

    def test_group_member_access(self, group):
        """Test that Group can access members when the Group is read-only.
        """
        assert len(group.members) == 2

    def test_group_moved_member_access(self, group, tmpdir):
        """Test that Group can access members when the Group is read-only,
        and when a member has been moved since the Group was last used with
        write-permissions.
        """
        with tmpdir.as_cwd():
            t = dtr.Treant('lark')
            t.location = 'somewhere/else'

        assert len(group.members) == 2
