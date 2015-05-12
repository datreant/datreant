"""Interface tests for Containers.

"""

import mdsynthesis as mds
import pandas as pd
import numpy as np
import pytest
import os
import shutil

containername = 'testcontainer'
simname = 'testsim'
groupname = 'testgroup'


class TestContainer:
    """Test generic Container features.

    """
    @pytest.fixture
    def container(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.containers.Container(containername)
        return c

    def test_init(self, container, tmpdir):
        """Test basic Container init.

        """
        assert container.name == containername
        assert container.location == tmpdir.strpath
        assert container.containertype == 'Container'
        assert container.basedir == os.path.join(tmpdir.strpath, containername)

    class TestTags:
        """Test container tags.

        """

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
        """Test container categories.

        """

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

    class TestData:
        """Test data storage and retrieval.

        """
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

            """Mixin class for pandas tests

            """
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

            """Test numpy datastructure storage and retrieval.

            """
            datafile = mds.core.persistence.npdatafile

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

            """Test pandas datastructure storage and retrieval.

            """
            datafile = mds.core.persistence.pydatafile

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


class TestSim:
    """Test Sim-specific features.

    """
    @pytest.fixture
    def sim(self, tmpdir):
        with tmpdir.as_cwd():
            s = mds.Sim(simname)
        return s

    def test_init(self, sim, tmpdir):
        """Test basic Container init.

        """
        assert sim.name == simname
        assert sim.location == tmpdir.strpath
        assert sim.containertype == 'Sim'

    def test_add_universe(self, sim):
        pass
