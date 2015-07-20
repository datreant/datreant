"""Tests for filesystem elements, including Foxhound and related objects.

"""

import datreant as dtr
import pytest
import os
import py.path


class TestFoxhound:
    """Test Foxhound functionality"""

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.treants.Container('testtreant')
        return c

    @pytest.fixture
    def sim(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.Treant('testsim')
        return c

    @pytest.fixture
    def group(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.Group('testgroup')
        return c
