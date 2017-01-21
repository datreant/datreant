"""Tests for filesystem elements, including Foxhound and related objects.

"""

import pytest

import datreant.core as dtr


class TestFoxhound:
    """Test Foxhound functionality"""

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.Treant('testtreant')
        return t

    @pytest.fixture
    def group(self, tmpdir):
        with tmpdir.as_cwd():
            g = dtr.Group('testgroup')
        return g
