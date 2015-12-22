"""Interface tests for Treants.

"""

import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py

from pytest import mark

import datreant as dtr
from datreant.backends.pytables import TreantFileHDF5
from datreant.backends.pyjson import TreantFileJSON
from datreant.backends.pyyaml import TreantFileYAML

def add_tags(treantfile, tags):
    treantfile.add_tags(*tags)

def add_tags_individually(treantfile, tags):
    for tag in tags:
        treantfile.add_tags(tag)

class TreantFileMixin:
    """Test generic TreantFile features"""

    def test_add_tags(self, treantfile):
        treantfile.add_tags('marklar')
        assert 'marklar' in treantfile.get_tags()

        treantfile.add_tags('lark', 'bark')
        assert 'marklar' in treantfile.get_tags()
        assert 'lark' in treantfile.get_tags()
        assert 'bark' in treantfile.get_tags()

#    def test_remove_tags(self, treantfile):
#        treantfile.tags.add('marklar')
#        assert 'marklar' in treantfile.tags
#        treantfile.tags.remove('marklar')
#        assert 'marklar' not in treantfile.tags
#
#        treantfile.tags.add('marklar')
#        treantfile.tags.add('lark', 'bark')
#        treantfile.tags.add(['fark', 'bark'])
#        assert 'marklar' in treantfile.tags
#        assert 'lark' in treantfile.tags
#        assert 'bark' in treantfile.tags
#        assert 'fark' in treantfile.tags
#        assert len(treantfile.tags) == 4
#
#        treantfile.tags.remove('fark')
#        assert 'fark' not in treantfile.tags
#        assert len(treantfile.tags) == 3
#        treantfile.tags.remove('fark')
#        assert len(treantfile.tags) == 3
#
#        treantfile.tags.remove(all=True)
#        assert len(treantfile.tags) == 0

    @mark.bench('add_tags')
    def test_add_single_tag(self, treantfile):
        tags = ['yes this {}'.format(i) for i in range(100)]
        treantfile.add_tags(*tags)

        add_tags(treantfile, 'a single tag')

    @mark.bench('add_tags')
    def test_add_many_tags_at_once(self, treantfile):
        tags = ['yes this {}'.format(i) for i in range(100)]
        add_tags(treantfile, tags)

    @mark.bench('add_tags_individually')
    def test_add_many_tags_individually(self, treantfile):
        tags = ['yes this {}'.format(i) for i in range(10)]
        add_tags_individually(treantfile, tags)
#
#    def test_add_categories(self, treantfile):
#        treantfile.categories.add(marklar=42)
#        assert 'marklar' in treantfile.categories
#
#        treantfile.categories.add({'bark': 'snark'}, lark=27)
#        assert 'bark' in treantfile.categories
#        assert 'snark' not in treantfile.categories
#        assert 'bark' in treantfile.categories
#
#        assert treantfile.categories['bark'] == 'snark'
#        assert treantfile.categories['lark'] == '27'
#
#        treantfile.categories['lark'] = 42
#        assert treantfile.categories['lark'] == '42'
#
#    def test_remove_categories(self, treantfile):
#        treantfile.categories.add(marklar=42)
#        assert 'marklar' in treantfile.categories
#
#        treantfile.categories.remove('marklar')
#        assert 'marklar' not in treantfile.categories
#
#        treantfile.categories.add({'bark': 'snark'}, lark=27)
#        del treantfile.categories['bark']
#        assert 'bark' not in treantfile.categories
#
#        # should just work, even if key isn't present
#        treantfile.categories.remove('smark')
#
#        treantfile.categories['lark'] = 42
#        treantfile.categories['fark'] = 32.3
#
#        treantfile.categories.remove(all=True)
#        assert len(treantfile.categories) == 0
#
#    def test_add_wrong(self, treantfile):
#        with pytest.raises(TypeError):
#            treantfile.categories.add('temperature', 300)
#
#        with pytest.raises(TypeError):
#            treantfile.categories.add(['mark', 'matt'])
#
#    def test_KeyError(self, treantfile):
#        with pytest.raises(KeyError):
#            treantfile.categories['hello?']


class TestTreantFileHDF5(TreantFileMixin):

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            c = TreantFileHDF5('testfile.h5')
        return c


class TestTreantFileJSON(TreantFileMixin):

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            c = TreantFileJSON('testfile.json')
        return c


class TestTreantFileYAML(TreantFileMixin):

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            c = TreantFileYAML('testfile.yaml')
        return c
