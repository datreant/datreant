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

        @pytest.fixture
        def collection(self, collection, tmpdir):
            with tmpdir.as_cwd():
                s1 = mds.Sim('lark')
                s2 = mds.Sim('hark')
                g3 = mds.Group('linus')

            collection.add(s1, [g3, s2])
            return collection

        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.uuid] = datastruct

                np.testing.assert_equal(collection.data.retrieve(self.handle),
                                        agg)
                np.testing.assert_equal(collection.data[self.handle],
                                        agg)

        class PanelMixin:
            """Mixin class for pandas structures that don't support
            MultiIndexing.

            """
            handle = 'testdata'

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.uuid] = datastruct

                stored = collection.data.retrieve(self.handle)
                for item in agg:
                    np.testing.assert_equal(stored[item].values,
                                            agg[item].values)

                stored = collection.data[self.handle]
                for item in agg:
                    np.testing.assert_equal(stored[item].values,
                                            agg[item].values)

        class MultiIndexMixin(DataMixin):
            """Mixin class for pandas structures that support MultiIndexes"""

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.uuid] = datastruct

                def dict2multiindex(agg):
                    agg_mi = None
                    for member in agg:
                        d = agg[member].copy(deep=True)
                        label = len(d.index)*[member]
                        index = pd.MultiIndex.from_arrays([label, d.index])
                        d.index = index

                        if agg_mi is not None:
                            agg_mi = agg_mi.append(d)
                        else:
                            agg_mi = d

                    return agg_mi

                np.testing.assert_equal(
                        collection.data.retrieve(self.handle).values,
                        dict2multiindex(agg).values)
                np.testing.assert_equal(
                        collection.data[self.handle].values,
                        dict2multiindex(agg).values)

        class Test_Series(test_data.Series, MultiIndexMixin):
            pass

        class Test_DataFrame(test_data.DataFrame, MultiIndexMixin):
            pass

        class Test_Blank_DataFrame(test_data.Blank_DataFrame, MultiIndexMixin):
            pass

        class Test_Wide_Blank_DataFrame(test_data.Wide_Blank_DataFrame,
                                        MultiIndexMixin):
            pass

        class Test_Thin_Blank_DataFrame(test_data.Thin_Blank_DataFrame,
                                        MultiIndexMixin):
            pass

        class Test_Panel(test_data.Panel, PanelMixin):
            pass

        class Test_Panel4D(test_data.Panel4D, PanelMixin):
            pass

        class Test_NumpyScalar(test_data.NumpyScalar, DataMixin):
            pass

        class Test_Numpy1D(test_data.Numpy1D, DataMixin):
            pass

        class Test_Numpy2D(test_data.Numpy2D, DataMixin):
            pass

        class Test_Wide_Numpy2D(test_data.Wide_Numpy2D, DataMixin):
            pass

        class Test_Thin_Numpy2D(test_data.Thin_Numpy2D, DataMixin):
            pass

        class Test_Numpy3D(test_data.Numpy3D, DataMixin):
            pass

        class Test_Numpy4D(test_data.Numpy4D, DataMixin):
            pass

        class Test_List(test_data.List, DataMixin):
            pass

        class Test_Dict(test_data.Dict, DataMixin):
            pass

        class Test_Tuple(test_data.Tuple, DataMixin):
            pass

        class Test_Set(test_data.Set, DataMixin):
            pass

        class Test_Dict_Mix(test_data.Dict_Mix, DataMixin):
            pass


class TestBundle(CollectionTests):
    """Test Bundle features"""

    @pytest.fixture
    def collection(self):
        return mds.Bundle()
