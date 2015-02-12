"""Interface tests for Containers. 

"""

import mdsynthesis as mds
import pytest
import os, shutil

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
        assert container.location == os.path.join(tmpdir.strpath, containername)
    
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

