"""Stress tests for file locks.

"""

import string
import multiprocessing as mp
import time
import pytest

from datreant.core import Treant


def pokefile(treantfilepath, string):
    """Add a number of tags to a Treant."""
    treant = Treant(treantfilepath)
    treant.tags.add(*["{}_{}".format(string, i) for i in range(100)])


def init_treant(tmpdir, tags):
    with tmpdir.as_cwd():
        tf = Treant('sprout', tags=tags)
    return tf


class TestTreantFile:

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = Treant('sprout')
        return t

    def test_death_by_1000_pokes(self, treant):
        pool = mp.Pool(processes=4)
        for i in range(10):
            pool.apply_async(pokefile, args=(treant.filepath,
                                             "run_{}".format(i)))
        pool.close()
        pool.join()

        assert len(treant.tags) == 1000

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
            tf = Treant('sprout')

        assert len(tf.tags) == num + 1
