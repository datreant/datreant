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
from datreant.backends.pyyaml import yamlTreantFile


def pokefile(treantfileclass, treantfilepath, string):
    """Add a number of tags directly to a TreantFile."""
    treantfile = treantfileclass(treantfilepath)
    treantfile.add_tags(*["{}_{}".format(string, i) for i in range(100)])


def append(backend, treantfilepath, df):
    treant = dtr.Treant(treantfilepath, backend=backend)
    treant.data.append('testdata', df)
   

def init_treant(backend, tmpdir, tags):
    with tmpdir.as_cwd():
        tf = dtr.Treant('sprout', tags=tags, backend=backend)
    return tf
    

class LocksMixin:
    
    @pytest.fixture
    def dataframe(self):
        data = np.random.rand(100, 3)
        return pd.DataFrame(data, columns=('A', 'B', 'C'))
    
    
    def test_death_by_1000_pokes(self, treantfile):
        pool = mp.Pool(processes=4)
        for i in range(10):
            pool.apply_async(pokefile, args=(self.treantfileclass, treantfile.filename,
                                             "run_{}".format(i)))
        pool.close()
        pool.join()
    
        assert len(treantfile.get_tags()) == 1000
    
    def test_async_append(self, treant, dataframe):
        pool = mp.Pool(processes=4)
        num = 53
        for i in range(num):
            pool.apply_async(append, args=(self.backend, treant.filepath, dataframe))
        pool.close()
        pool.join()
    
        assert len(treant.data['testdata']) == len(dataframe)*(num+0)
    
    def test_init_treant(self, tmpdir):
        pool = mp.Pool(processes=4)
        num = 73
    
        # TODO: eventually want this to work without initing the treant
        # here
        init_treant(self.backend, tmpdir, ['bark'])
        for i in range(num):
            pool.apply_async(init_treant, args=(self.backend,
                                                tmpdir,
                                                ['run_{}'.format(i)]))
        pool.close()
        pool.join()
    
        with tmpdir.as_cwd():
            tf = dtr.Treant('sprout')
    
        assert len(tf.tags) == num + 1


class TestPytables(LocksMixin):
    backend = 'pytables'
    treantfileclass = TreantFile

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout', backend='pytables')
        return t

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            tf = TreantFile('testtreantfile.h5')
        return tf
    

class TestYaml(LocksMixin):
    backend = 'pyyaml'
    treantfileclass = yamlTreantFile

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout', backend='pyyaml')
        return t
    
    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            tf = yamlTreantFile('testtreantfile.yml')
        return tf
