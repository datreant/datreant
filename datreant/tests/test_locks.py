"""Stress tests for file locks.

"""

import string
import multiprocessing as mp
import time
import pytest

import datreant as dtr

from datreant.backends.statefiles import TreantFile


def pokefile(treantfileclass, treantfilepath, string):
    """Add a number of tags directly to a TreantFile."""
    treantfile = treantfileclass(treantfilepath)
    treantfile.add_tags(*["{}_{}".format(string, i) for i in range(100)])


def init_treant(tmpdir, tags):
    with tmpdir.as_cwd():
        tf = dtr.Treant('sprout', tags=tags)
    return tf


class TestTreantFile:

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout')
        return t

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            tf = TreantFile('testtreantfile.json')
        return tf

    def test_death_by_1000_pokes(self, treantfile):
        pool = mp.Pool(processes=4)
        for i in range(10):
            pool.apply_async(pokefile, args=(TreantFile,
                                             treantfile.filename,
                                             "run_{}".format(i)))
        pool.close()
        pool.join()

        assert len(treantfile.get_tags()) == 1000

    def test_init_treant(self, tmpdir):
        pool = mp.Pool(processes=4)
        num = 73

        # TODO: eventually want this to work without initing the treant
        # here
        init_treant(tmpdir, ['bark'])
        for i in range(num):
            pool.apply_async(init_treant, args=(tmpdir,
                                                ['run_{}'.format(i)]))
        pool.close()
        pool.join()

        with tmpdir.as_cwd():
            tf = dtr.Treant('sprout')

        assert len(tf.tags) == num + 1
