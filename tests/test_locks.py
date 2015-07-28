"""Stress tests for file locks.

"""

import random
import string
import multiprocessing as mp
import time
import pytest

from datreant.persistence import TreantFile

@pytest.fixture
def treantfile(tmpdir):
    with tmpdir.as_cwd():
        tf = TreantFile('testtreantfile.h5')
    return tf

def pokefile(treantfilepath, string):
    """Add two tags, then delete a tag, directly from a TreantFile"""
    treantfile = TreantFile(treantfilepath)

    #for i in range(2):
    #    tag = "{}_{}".format(string, i)
    #    treantfile.add_tags(tag)

    treantfile._open_w()
    treantfile.logger.info("task {}, fd {}".format(string, treantfile.handle.fileno()))
    treantfile._close()

    #treantfile.add_tags("task_{}".format(string))

    treantfile.add_tags(*["{}_{}".format(string, i) for i in range(1000)])

    #treantfile.del_tags(random.choice(list(treantfile.get_tags())))

def test_death_by_1000_pokes(treantfile):
    #for i in xrange(2):
    #    pokefile(treantfile.filename, "run_{}".format(i))
    #mp.set_start_method('spawn')
    pool = mp.Pool(processes=4)
    for i in range(12):
        pool.apply_async(pokefile, args=(treantfile.filename, 
                                         "run_{}".format(i)))
    pool.close()
    pool.join()

    assert len(treantfile.get_tags()) == 12000
    #assert len(treantfile.get_tags()) == 1000
