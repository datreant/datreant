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
        class TestPandas:
            """Test pandas datastructure storage and retrieval.

            """
            @pytest.fixture
            def series(self):
                data = np.random.rand(10000)
                s = pd.Series(data)
                return s
    
            @pytest.fixture
            def dataframe(self):
                data = np.random.rand(10000,3)
                df = pd.DataFrame(data, columns=('A', 'B', 'C'))
                return df

            @pytest.fixture
            def panel(self):
                data = np.random.rand(4,10000,3)
                p = pd.Panel(data, items=('I', 'II', 'III', 'IV'), minor_axis=('A', 'B', 'C'))
                return p

            def test_lock_series(series)

        class TestNumpy:
            """Test pandas datastructure storage and retrieval.

            """
            def test_lock

        class TestPython:
            """Test pandas datastructure storage and retrieval.

            """
            def test_lock




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

