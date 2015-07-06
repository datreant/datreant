"""Tests for Bundle.

"""

import mdsynthesis as mds
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py
import test_data

import MDAnalysis
from MDAnalysisTests.datafiles import GRO, XTC

def do_stuff(cont):
    return cont.name + cont.uuid

def return_nothing(cont):
    b = cont.name + cont.uuid
    
class CollectionTests():
    """Tests for common elements of Group.members and Bundle"""

    def test_add_members(self, collection, tmpdir):
        """Try adding members in a number of ways"""
        with tmpdir.as_cwd():
            s1 = mds.Sim('lark')
            s2 = mds.Sim('hark')
            g3 = mds.Group('linus')
    
            collection.add(s1, [g3, s2])
    
            for cont in (s1, s2, g3):
                assert cont in collection
    
            s4 = mds.Sim('snoopy')
            collection.add([[s4], s2])
            assert s4 in collection

    def test_get_members(self, collection, tmpdir):
        """Access members with indexing and slicing"""
        with tmpdir.as_cwd():
            s1 = mds.Sim('larry')
            g2 = mds.Group('curly')
            s3 = mds.Sim('moe')
    
            collection.add([[[s1, [g2, [s3]]]]])
    
            assert collection[1] == g2
    
            c4 = mds.containers.Container('shemp')
            collection.add(c4)
    
            for member in (s1, g2, s3):
                assert member in collection[:3]
    
            assert c4 not in collection[:3]
            assert c4 == collection[-1]
    
    def test_remove_members(self, collection, tmpdir):
        """Try removing members"""
        with tmpdir.as_cwd():
            g1 = mds.Group('lion-o')
            s2 = mds.Sim('cheetara')
            s3 = mds.Container('snarf')
    
            collection.add(s3, g1, s2)
    
            for cont in (g1, s2, s3):
                assert cont in collection
    
            collection.remove(1)
            assert g1 not in collection
    
            collection.remove(s2)
            assert s2 not in collection
    
    def test_member_attributes(self, collection, tmpdir):
        """Get member uuids, names, and containertypes"""
        with tmpdir.as_cwd():
            c1 = mds.containers.Container('bigger')
            g2 = mds.Group('faster')
            s3 = mds.Sim('stronger')
    
        collection.add(c1, g2, s3)
    
        uuids = [cont.uuid for cont in [c1, g2, s3]]
        assert collection.uuids == uuids
    
        names = [cont.name for cont in [c1, g2, s3]]
        assert collection.names == names
    
        containertypes = [cont.containertype for cont in [c1, g2, s3]]
        assert collection.containertypes == containertypes
    
    def test_map(self, collection, tmpdir):
        with tmpdir.as_cwd():
            s1 = mds.Sim('lark')
            s2 = mds.Sim('hark')
            g3 = mds.Group('linus')

        collection.add(s1, s2, g3)
    
        comp = [cont.name + cont.uuid for cont in collection]
        assert collection.map(do_stuff) == comp
        assert collection.map(do_stuff, processes=2) == comp
    
        assert collection.map(return_nothing) is None
        assert collection.map(return_nothing, processes=2) is None
    
    class TestMemberData:
        """Test member data functionality"""
        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_retrieve_data(self, container, datastruct):
                container.data.add(self.handle, datastruct)
                np.testing.assert_equal(container.data.retrieve(self.handle),
                                        datastruct)
                np.testing.assert_equal(container.data[self.handle],
                                        datastruct)

        class PandasMixin(DataMixin):
            """Mixin class for pandas tests"""
            datafile = mds.core.persistence.pddatafile

            def test_retrieve_data(self, container, datastruct):
                container.data.add(self.handle, datastruct)
                np.testing.assert_equal(
                        container.data.retrieve(self.handle).values,
                        datastruct.values)
                np.testing.assert_equal(
                        container.data[self.handle].values,
                        datastruct.values)

        class Test_Series(test_data.Series, PandasMixin):
            pass

        class Test_DataFrame(test_data.DataFrame, PandasMixin):
            pass

        class Test_Blank_DataFrame(test_data.Blank_DataFrame, PandasMixin):
            pass

        class Test_Wide_Blank_DataFrame(test_data.Wide_Blank_DataFrame, PandasMixin):
            pass

        class Test_Thin_Blank_DataFrame(test_data.Thin_Blank_DataFrame, PandasMixin):
            pass

        class Test_Panel(test_data.Panel, PandasMixin):
            pass

        class Test_Panel4D(test_data.Panel4D, PandasMixin):
            pass

        class NumpyMixin(DataMixin):
            """Test numpy datastructure storage and retrieval"""
            datafile = mds.core.persistence.npdatafile

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
            datafile = mds.core.persistence.pydatafile

            def test_overwrite_data(self, container, datastruct):
                container.data[self.handle] = datastruct

                # overwrite the data with a scalar
                container.data[self.handle] = 23
                assert container.data[self.handle] == 23

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

class TestBundle:
    """Test Bundle features"""

    @pytest.fixture
    def sim_data(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.containers.Container(TestContainer.containername)
        return c

    @pytest.fixture
    def group_data(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.containers.Container(TestContainer.containername)
        return c

    def test_init(self, container, tmpdir):
        """Test basic Container init"""
        assert container.name == self.containername
        assert container.location == tmpdir.strpath
        assert container.containertype == self.containertype
        assert container.basedir == os.path.join(tmpdir.strpath,
                                                 self.containername)

    def test_cmp(self, tmpdir):
        """Test the comparison of Containers when sorting"""
        with tmpdir.as_cwd():
            c1 = self.containerclass('a')
            c2 = self.containerclass('b')
            c3 = self.containerclass('c')

        assert sorted([c3, c2, c1]) == [c1, c2, c3]
        assert c1 <= c2 < c3
        assert c3 >= c2 > c1

    class TestTags:
        """Test container tags"""

        def test_add_tags(self, container):
            container.tags.add('marklar')
            assert 'marklar' in container.tags

            container.tags.add('lark', 'bark')
            assert 'marklar' in container.tags
            assert 'lark' in container.tags
            assert 'bark' in container.tags

        def test_remove_tags(self, container):
            container.tags.add('marklar')
            assert 'marklar' in container.tags
            container.tags.remove('marklar')
            assert 'marklar' not in container.tags

            container.tags.add('marklar')
            container.tags.add('lark', 'bark')
            container.tags.add(['fark', 'bark'])
            assert 'marklar' in container.tags
            assert 'lark' in container.tags
            assert 'bark' in container.tags
            assert 'fark' in container.tags
            assert len(container.tags) == 4

            container.tags.remove('fark')
            assert 'fark' not in container.tags
            assert len(container.tags) == 3
            container.tags.remove('fark')
            assert len(container.tags) == 3

            container.tags.remove(all=True)
            assert len(container.tags) == 0

    class TestCategories:
        """Test container categories"""

        def test_add_categories(self, container):
            container.categories.add(marklar=42)
            assert 'marklar' in container.categories

            container.categories.add({'bark': 'snark'}, lark=27)
            assert 'bark' in container.categories
            assert 'snark' not in container.categories
            assert 'bark' in container.categories

            assert container.categories['bark'] == 'snark'
            assert container.categories['lark'] == '27'

            container.categories['lark'] = 42
            assert container.categories['lark'] == '42'

        def test_remove_categories(self, container):
            container.categories.add(marklar=42)
            assert 'marklar' in container.categories

            container.categories.remove('marklar')
            assert 'marklar' not in container.categories

            container.categories.add({'bark': 'snark'}, lark=27)
            del container.categories['bark']
            assert 'bark' not in container.categories

            container.categories['lark'] = 42
            container.categories['fark'] = 32.3

            container.categories.remove(all=True)
            assert len(container.categories) == 0

        def test_add_wrong(self, container):
            with pytest.raises(TypeError):
                container.categories.add('temperature', 300)

            with pytest.raises(TypeError):
                container.categories.add(['mark', 'matt'])

        def test_KeyError(self, container):
            with pytest.raises(KeyError):
                container.categories['hello?']

    class TestData:
        """Test data storage and retrieval"""

        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_add_data(self, container, datastruct):
                container.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(container.basedir,
                                                   self.handle,
                                                   self.datafile))

            def test_remove_data(self, container, datastruct):
                container.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(container.basedir,
                                                   self.handle,
                                                   self.datafile))

                container.data.remove('testdata')
                assert not os.path.exists(os.path.join(container.basedir,
                                                       self.handle,
                                                       self.datafile))

        class PandasMixin(DataMixin):
            """Mixin class for pandas tests"""
            datafile = mds.core.persistence.pddatafile

        class Test_Series(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                data = np.random.rand(10000)
                return pd.Series(data)

        class Test_DataFrame(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                data = np.random.rand(10000, 3)
                return pd.DataFrame(data, columns=('A', 'B', 'C'))

        class Test_Blank_DataFrame(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                return pd.DataFrame(np.zeros((10, 10)))

        class Test_Wide_Blank_DataFrame(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                return pd.DataFrame(np.zeros((1, 10)))

        class Test_Thin_Blank_DataFrame(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                return pd.DataFrame(np.zeros((10, 1)))

        class Test_Panel(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                data = np.random.rand(4, 10000, 3)
                return pd.Panel(data, items=('I', 'II', 'III', 'IV'),
                                minor_axis=('A', 'B', 'C'))

        class Test_Panel4D(PandasMixin):
            @pytest.fixture
            def datastruct(self):
                data = np.random.rand(2, 4, 10000, 3)
                return pd.Panel4D(data, labels=('gallahad', 'lancelot'),
                                  items=('I', 'II', 'III', 'IV'),
                                  minor_axis=('A', 'B', 'C'))

        class NumpyMixin(DataMixin):
            """Test numpy datastructure storage and retrieval"""
            datafile = mds.core.persistence.npdatafile

        class Test_NumpyScalar(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.array(20)

        class Test_Numpy1D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.random.rand(10000)

        class Test_Numpy2D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.random.rand(10000, 500)

        class Test_Wide_Numpy2D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.zeros((1, 10))

        class Test_Thin_Numpy2D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.zeros((10, 1))

        class Test_Numpy3D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.random.rand(4, 10000, 45)

        class Test_Numpy4D(NumpyMixin):
            @pytest.fixture
            def datastruct(self):
                return np.random.rand(2, 4, 10000, 45)

        class PythonMixin(DataMixin):
            """Test pandas datastructure storage and retrieval"""
            datafile = mds.core.persistence.pydatafile

            def test_overwrite_data(self, container, datastruct):
                container.data[self.handle] = datastruct

                # overwrite the data with a scalar
                container.data[self.handle] = 23
                assert container.data[self.handle] == 23

        class Test_List(PythonMixin):
            @pytest.fixture
            def datastruct(self):
                return ['galahad', 'lancelot', 42, 'arthur', 3.14159]

        class Test_Dict(PythonMixin):
            @pytest.fixture
            def datastruct(self):
                return {'pure': 'galahad', 'brave': 'lancelot', 'answer': 42,
                        'king': 'arthur', 'pi-ish': 3.14159}

        class Test_Tuple(PythonMixin):
            @pytest.fixture
            def datastruct(self):
                return ('arthur', 3.14159)

        class Test_Set(PythonMixin):
            @pytest.fixture
            def datastruct(self):
                return {'arthur', 3.14159, 'seahorses'}

        class Test_Dict_Mix(PythonMixin):
            @pytest.fixture
            def datastruct(self):
                return {'an array': np.random.rand(100, 46),
                        'another': np.random.rand(3, 45, 2)}

