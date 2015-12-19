"""Stress tests for file locks.

"""

import random
import string
import multiprocessing as mp
import time
import pytest
import numpy as np
import pandas as pd

import datreant as dtr

from datreant.backends.pytables import TreantFile


class Locks:
    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout')
        return t
    
    @pytest.fixture
    def dataframe(self):
        data = np.random.rand(100, 3)
        return pd.DataFrame(data, columns=('A', 'B', 'C'))
    
    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            tf = TreantFile('testtreantfile.h5')
        return tf
    
    def pokefile(self, treantfilepath, string):
        """Add a number of tags directly to a TreantFile."""
        treantfile = TreantFile(treantfilepath)
        treantfile.add_tags(*["{}_{}".format(string, i) for i in range(1000)])
    
    def test_death_by_99000_pokes(self, treantfile):
        pool = mp.Pool(processes=4)
        for i in range(99):
            pool.apply_async(self.pokefile, args=(treantfile.filename,
                                             "run_{}".format(i)))
        pool.close()
        pool.join()
    
        assert len(treantfile.get_tags()) == 99000
    
    def append(self, treantfilepath, df):
        treant = dtr.Treant(treantfilepath)
        treant.data.append('testdata', df)
    
    def test_async_append(self, treant, dataframe):
        pool = mp.Pool(processes=4)
        num = 53
        for i in range(num):
            pool.apply_async(self.append, args=(treant.filepath, dataframe))
        pool.close()
        pool.join()
    
        assert len(treant.data['testdata']) == len(dataframe)*(num+0)
    
    def init_treant(self, tmpdir, tags):
        with tmpdir.as_cwd():
            tf = dtr.Treant('sprout', tags=tags)
        return tf
    
    def test_init_treant(self, tmpdir):
        pool = mp.Pool(processes=4)
        num = 73
    
        # TODO: eventually want this to work without initing the treant
        # here
        self.init_treant(tmpdir, ['bark'])
        for i in range(num):
            pool.apply_async(self.init_treant, args=(tmpdir, ['run_{}'.format(i)]))
        pool.close()
        pool.join()
    
        with tmpdir.as_cwd():
            tf = dtr.Treant('sprout')
    
        assert len(tf.tags) == num + 1


class TestPytables(Locks):
    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout', backend='pytables')
        return t
    

class TestYaml(Locks):
    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout', backend='yaml')
        return t
    
