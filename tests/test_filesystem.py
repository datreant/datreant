"""Tests for filesystem elements, including Foxhound and related objects.

"""

import mdsynthesis as mds
import pytest
import os
import py.path

import MDAnalysis
from MDAnalysisTests.datafiles import GRO, XTC


class TestFoxhound:
    """Test Foxhound functionality"""

    @pytest.fixture
    def container(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.containers.Container('testcontainer')
        return c

    @pytest.fixture
    def sim(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.Sim('testsim')
        return c

    @pytest.fixture
    def group(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.Group('testgroup')
        return c


class TestUniversehound:
    """Test Universehound functionality"""
    GRO = 'md.gro'
    XTC = 'md.xtc'

    @pytest.fixture
    def sim_external(self, tmpdir):
        """Create a Sim object with an externally-stored set of
        universe files.
        """
        with tmpdir.as_cwd():
            s = mds.Sim('testsim')
            sub = tmpdir.mkdir('sub')
            GRO_t = sub.join(self.GRO)
            XTC_t = sub.join(self.XTC)
            py.path.local(GRO).copy(GRO_t)
            py.path.local(XTC).copy(XTC_t)

            s.universes.add('spam', GRO_t.strpath, XTC_t.strpath)
        return s

    @pytest.fixture
    def sim_internal(self, tmpdir):
        """Create a Sim object with an internally-stored set of
        universe files.
        """
        with tmpdir.as_cwd():
            s = mds.Sim('testsim')
            sub = py.path.local(s.basedir).mkdir('sub')
            GRO_t = sub.join(self.GRO)
            XTC_t = sub.join(self.XTC)
            py.path.local(GRO).copy(GRO_t)
            py.path.local(XTC).copy(XTC_t)

            s.universes.add('spam', GRO_t.strpath, XTC_t.strpath)
        return s

    def test_move_Sim_external(self, sim_external, tmpdir):
        """Test that when we move the Sim object, it can still
        find its externally-stored set of universe files.
        """
        sim_external.location = tmpdir.mkdir('test').strpath
        assert isinstance(sim_external.universe, MDAnalysis.Universe)

        top, traj = sim_external.universes.define('spam')
        assert sim_external.universe.filename == top
        assert sim_external.universe.trajectory.filename == traj[0]

    def test_move_Sim_internal(self, sim_internal, tmpdir):
        sim_internal.location = tmpdir.mkdir('test').strpath
        assert isinstance(sim_internal.universe, MDAnalysis.Universe)

        top, traj = sim_internal.universes.define('spam')
        assert sim_internal.universe.filename == top
        assert sim_internal.universe.trajectory.filename == traj[0]

    def test_move_Sim_internal_copy(self, sim_internal, tmpdir):
        """We move the whole Sim, including its internal trajectory files, and
        we put a copy of one trajectory file where the internal ones used to
        be.
        """
        grofile = tmpdir.join(sim_internal.name, 'sub', self.GRO)
        assert sim_internal.universes.define('spam')[0] == grofile.strpath

        sim_internal.location = tmpdir.mkdir('test').strpath
        tmpdir.join('test', sim_internal.name, 'sub', self.GRO).copy(
                tmpdir.mkdir(sim_internal.name).mkdir('sub'))

        with pytest.raises(IOError):
            sim_internal.universe
