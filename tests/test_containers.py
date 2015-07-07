"""Interface tests for Containers.

"""

import mdsynthesis as mds
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py
import test_data
import test_bundle

import MDAnalysis
from MDAnalysisTests.datafiles import GRO, XTC


class TestContainer:
    """Test generic Container features"""
    containername = 'testcontainer'
    containertype = 'Container'
    containerclass = mds.containers.Container

    @pytest.fixture
    def container(self, tmpdir):
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

    def test_regen(self, tmpdir):
        """Test regenerating Container.

        - create Container
        - modify Container a little
        - create same Container (should regenerate)
        - check that modifications were saved
        """
        with tmpdir.as_cwd():
            C1 = self.containerclass('regen')
            C1.tags.add('fantastic')
            C2 = self.containerclass('regen')  # should be regen of C1
            assert 'fantastic' in C2.tags

            # they point to the same file, but they are not the same object
            assert C1 is not C2

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


class TestSim(TestContainer):
    """Test Sim-specific features"""
    containername = 'testsim'
    containertype = 'Sim'
    containerclass = mds.Sim

    @pytest.fixture
    def container(self, tmpdir):
        with tmpdir.as_cwd():
            s = mds.Sim(TestSim.containername)
        return s

    class TestUniverses:
        """Test universe functionality"""

        def test_add_universe(self, container):
            """Test adding new unverse definitions"""
            container.universes.add('spam', GRO, XTC)

            assert 'spam' in container.universes
            assert isinstance(container.universes['spam'], MDAnalysis.Universe)

            assert container.universe.filename == GRO
            assert container.universe.trajectory.filename == XTC

        def test_remove_universe(self, container):
            """Test universe removal"""
            container.universes.add('spam', GRO, XTC)
            container.universes.add('eggs', GRO, XTC)

            assert 'spam' in container.universes

            container.universes.remove('spam')

            assert 'spam' not in container.universes
            assert 'eggs' in container.universes

            with pytest.raises(KeyError):
                container.universes.remove('ham')

        def test_set_default_universe(self, container):
            """Test that a default universe exists, and that it's settable"""
            container.universes.add('lolcats', GRO, XTC)

            assert container.universe.filename == GRO
            assert container.universe.trajectory.filename == XTC
            assert container._uname == 'lolcats'
            container.universes.deactivate()

            container.universes.add('megaman', GRO)
            container.universes.default('megaman')

            assert container.universes.default() == 'megaman'

            assert container.universe.filename == GRO
            assert container.universe.trajectory.filename == GRO
            assert container._uname == 'megaman'

            container.universes.remove('megaman')

            assert container.universes.default() == None

        def test_set_resnums(self, container):
            """Test that we can add resnums to a universe."""
            container.universes.add('lolcats', GRO, XTC)

            protein = container.universe.selectAtoms('protein')
            resids = protein.residues.resids()
            protein.residues.set_resnum(resids + 3)

            container.universes.resnums('lolcats',
                                        container.universe.atoms.resnums())

            container.universes['lolcats']

            protein = container.universe.selectAtoms('protein')
            assert (resids + 3 == protein.residues.resnums()).all()

            # BUG IN MDANALYSIS PREVENTS RESETTING OF RESNUMS
            # protein.residues.set_resnum(resids + 6)

            # assert (protein.residues.resnums == resids + 6).all()
            # container.universes.resnums('lolcats',
            #                            container.universe.atoms.resnums())

            # container.universes['lolcats']

            # protein = container.universe.selectAtoms('protein')
            # assert (resids + 6 == protein.residues.resnums()).all()

        def test_KeyError(self, container):
            """Test that a KeyError raised when trying to activate a Universe
            that doesn't exist.
            """
            with pytest.raises(KeyError):
                container.universes.activate('ham')

            container.universes.add('lolcats', GRO, XTC)

            with pytest.raises(KeyError):
                container.universes.activate('eggs')

    class TestSelections:
        """Test stored selections functionality"""
        @pytest.fixture
        def container(self, tmpdir):
            with tmpdir.as_cwd():
                s = mds.Sim(TestSim.containername)
                s.universes.add('spam', GRO, XTC)
            return s

        def test_add_selection(self, container):
            """Test adding new selection definitions"""

            container.selections.add('CA', 'protein and name CA')
            container.selections.add('someres', 'resid 12')

            CA = container.universe.selectAtoms('protein and name CA')
            someres = container.universe.selectAtoms('resid 12')

            assert (CA.indices() == container.selections['CA'].indices()).all()
            assert (someres.indices() ==
                    container.selections['someres'].indices()).all()

        def test_remove_selection(self, container):
            """Test universe removal"""

            container.selections.add('CA', 'protein and name CA')
            container.selections.add('someres', 'resid 12')

            assert 'CA' in container.selections
            assert 'someres' in container.selections

            container.selections.remove('CA')

            assert 'CA' not in container.selections
            assert 'someres' in container.selections

            del container.selections['someres']
            container.selections.add('moreres', 'resid 12:20')

            assert 'CA' not in container.selections
            assert 'someres' not in container.selections
            assert 'moreres' in container.selections

            with pytest.raises(KeyError):
                container.selections.remove('someres')

            with pytest.raises(KeyError):
                del container.selections['CA']

        def test_selection_keys(self, container):
            container.selections.add('CA', 'protein and name CA')
            container.selections.add('someres', 'resid 12')

            assert set(('CA', 'someres')) == set(container.selections.keys())

        def test_selection_define(self, container):
            CA = 'protein and name CA'
            container.selections.add('CA', CA)

            assert container.selections.define('CA')[0] == CA

        def test_selection_get(self, container):
            with pytest.raises(KeyError):
                container.selections['CA']

        def test_add_selections_multiple_strings_via_add(self, container):
            """Add a selection that has multiple selection strings"""
            container.selections.add('funky town', 'name N', 'name CA')
            assert 'funky town' in container.selections

            ref = container.universe.selectAtoms('name N', 'name CA')
            sel = container.selections['funky town']
            assert (ref.indices() == sel.indices()).all()

        def test_add_selections_multiple_strings_via_setitem(self, container):
            """Add a selection that has multiple selection strings"""
            container.selections['funky town 2'] = 'name N', 'name CA'
            assert 'funky town 2' in container.selections

            ref = container.universe.selectAtoms('name N', 'name CA')
            sel = container.selections['funky town 2']
            assert (ref.indices() == sel.indices()).all()

        def test_add_selection_as_atomgroup_via_add(self, container):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = container.universe.atoms[:10:2]

            container.selections.add('ag sel', ag)
            assert 'ag sel' in container.selections

            ag2 = container.selections['ag sel']
            assert (ag.indices() == ag2.indices()).all()

        def test_add_selection_as_atomgroup_via_setitem(self, container):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = container.universe.atoms[25:50:3]

            container.selections['ag sel 2'] = ag
            assert 'ag sel 2' in container.selections

            ag2 = container.selections['ag sel 2']
            assert (ag.indices() == ag2.indices()).all()


class TestGroup(TestContainer):
    """Test Group-specific features.

    """
    containername = 'testgroup'
    containertype = 'Group'
    containerclass = mds.Group

    @pytest.fixture
    def container(self, tmpdir):
        with tmpdir.as_cwd():
            g = mds.Group(TestGroup.containername)
        return g

    def test_repr(self, container):
        pass

    class TestMembers(test_bundle.CollectionTests):
        """Test member functionality"""

        @pytest.fixture
        def collection(self, container):
            return container.members


class TestReadOnly:
    """Test Container functionality when read-only"""

    GRO = 'md.gro'
    XTC = 'md.xtc'

    @pytest.fixture
    def container(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = mds.containers.Container('testcontainer')
            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    @pytest.fixture
    def sim(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = mds.Sim('testsim')

            # copy universe files to within the Sim's tree
            sub = py.path.local(c.basedir).mkdir('sub')
            GRO_t = sub.join(self.GRO)
            XTC_t = sub.join(self.XTC)
            py.path.local(GRO).copy(GRO_t)
            py.path.local(XTC).copy(XTC_t)

            c.universes.add('main', GRO_t.strpath, XTC_t.strpath)

            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    @pytest.fixture
    def group(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = mds.Group('testgroup')
            c.members.add(mds.Sim('lark'), mds.Group('bark'))
            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    def test_sim_universe_access(self, sim):
        """Test that Sim can access Universe when read-only.
        """
        assert isinstance(sim.universe, MDAnalysis.Universe)

    def test_sim_moved_universe_access(self, sim, tmpdir):
        """Test that Sim can access Universe when read-only, especially
        when universe files have moved with it (stale paths).
        """
        py.path.local(sim.basedir).chmod(0770, rec=True)
        sim.location = tmpdir.mkdir('test').strpath
        py.path.local(sim.basedir).chmod(0550, rec=True)

        assert isinstance(sim.universe, MDAnalysis.Universe)

        py.path.local(sim.basedir).chmod(0770, rec=True)

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
            s = mds.Sim('lark')
            s.location = 'somewhere/else'

        assert len(group.members) == 2
