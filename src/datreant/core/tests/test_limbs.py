import datreant.core as dtr

import pytest


class TestCategory(object):
    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            yield dtr.Treant('this')

    def test_setting_to_None_VE(self, treant):
        with pytest.raises(ValueError) as err:
            treant.categories['colour'] = None
        assert "Cannot set to 'None'" in str(err)

    def test_setting_None_to_delete(self, treant):
        treant.categories['colour'] = 'blue'
        treant.categories['size'] = 'large'

        treant.categories['size'] = None

        assert 'size' not in treant.categories
        assert 'colour' in treant.categories
