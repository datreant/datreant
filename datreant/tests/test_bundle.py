"""Tests for Bundle.

"""

import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py

import datreant as dtr
from datreant.tests import test_data


def do_stuff(cont):
    return cont.name + cont.uuid


def return_nothing(cont):
    b = cont.name + cont.uuid


class CollectionTests():
    """Tests for common elements of Group.members and Bundle"""

    @pytest.fixture
    def testtreant(self, tmpdir, request):
        with tmpdir.as_cwd():
            t = dtr.Treant('dummytreant')
        return t

    @pytest.fixture
    def testgroup(self, tmpdir, request):
        with tmpdir.as_cwd():
            g = dtr.Group('dummygroup')
            g.members.add(dtr.Treant('bark'), dtr.Treant('leaf'))
        return g

    def test_additive(self, tmpdir, testtreant, testgroup, collection):
        """Test that addition of treants and collections give Bundles.

        """
        with tmpdir.as_cwd():
            assert isinstance(testtreant + testgroup, dtr.Bundle)
            assert len(testtreant + testgroup) == 2

            # subtle, but important; Group.members is a collection,
            # while Group is a treant
            assert len(testtreant + testgroup.members) != 2
            assert (len(testtreant + testgroup.members) ==
                    len(testgroup.members) + 1)

            assert isinstance(testtreant + testgroup.members, dtr.Bundle)

            b = collection + testtreant + testgroup

            # beating a dead horse
            assert len(b) == 2
            assert (len(b + testgroup.members) ==
                    len(b) + len(testgroup.members))
            assert isinstance(b + testgroup.members, dtr.Bundle)

    def test_add_members(self, collection, tmpdir):
        """Try adding members in a number of ways"""
        with tmpdir.as_cwd():
            s1 = dtr.Treant('lark')
            s2 = dtr.Treant('hark')
            g3 = dtr.Group('linus')

            collection.add(s1, [g3, s2])

            for cont in (s1, s2, g3):
                assert cont in collection

            s4 = dtr.Treant('snoopy')
            collection.add([[s4], s2])
            assert s4 in collection

            # the group won't add members it alrady has
            # (operates as an ordered set)
            assert len(collection) == 4

    def test_add_members_glob(self, collection, tmpdir):
        """Try adding members with globbing"""
        with tmpdir.as_cwd():
            t1 = dtr.Treant('lark')
            t2 = dtr.Treant('hark')
            g3 = dtr.Group('linus')

            collection.add('*ark')

            for treant in (t1, t2):
                assert treant in collection

            assert g3 not in collection

    def test_get_members(self, collection, tmpdir):
        """Access members with indexing and slicing"""
        with tmpdir.as_cwd():
            s1 = dtr.Treant('larry')
            g2 = dtr.Group('curly')
            s3 = dtr.Treant('moe')

            collection.add([[[s1, [g2, [s3]]]]])

            assert collection[1] == g2

            c4 = dtr.treants.Treant('shemp')
            collection.add(c4)

            for member in (s1, g2, s3):
                assert member in collection[:3]

            assert c4 not in collection[:3]
            assert c4 == collection[-1]

    def test_remove_members(self, collection, tmpdir):
        """Try removing members"""
        with tmpdir.as_cwd():
            g1 = dtr.Group('lion-o')
            s2 = dtr.Treant('cheetara')
            s3 = dtr.Treant('snarf')

            collection.add(s3, g1, s2)

            for cont in (g1, s2, s3):
                assert cont in collection

            collection.remove(1)
            assert g1 not in collection

            collection.remove(s2)
            assert s2 not in collection

    def test_remove_members_name(self, collection, tmpdir):
        """Try removing members with names and globbing"""
        with tmpdir.as_cwd():
            t1 = dtr.Treant('lark')
            t2 = dtr.Treant('elsewhere/lark')
            t3 = dtr.Treant('hark')
            g = dtr.Group('linus')

            stuff = [t1, t2, t3, g]

            # test removal by name
            collection.add(stuff)
            for item in stuff:
                assert item in collection

            # should remove both treants with name 'lark'
            collection.remove('lark')

            for item in (t3, g):
                assert item in collection

            for item in (t1, t2):
                assert item not in collection

            # test removal by a unix-style glob pattern
            collection.add(stuff)
            for item in stuff:
                assert item in collection

            # should remove 'lark' and 'hark' treants
            collection.remove('*ark')

            assert g in collection

            for item in (t1, t2, t3):
                assert item not in collection

    def test_member_attributes(self, collection, tmpdir):
        """Get member uuids, names, and treanttypes"""
        with tmpdir.as_cwd():
            c1 = dtr.treants.Treant('bigger')
            g2 = dtr.Group('faster')
            s3 = dtr.Treant('stronger')

        collection.add(c1, g2, s3)

        uuids = [cont.uuid for cont in [c1, g2, s3]]
        assert collection.uuids == uuids

        names = [cont.name for cont in [c1, g2, s3]]
        assert collection.names == names

        treanttypes = [cont.treanttype for cont in [c1, g2, s3]]
        assert collection.treanttypes == treanttypes

    def test_map(self, collection, tmpdir):
        with tmpdir.as_cwd():
            s1 = dtr.Treant('lark')
            s2 = dtr.Treant('hark')
            g3 = dtr.Group('linus')

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
                s1 = dtr.Treant('lark')
                s2 = dtr.Treant('hark')
                g3 = dtr.Group('linus')

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
        return dtr.Bundle()
