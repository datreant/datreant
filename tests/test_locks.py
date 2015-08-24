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

    treantfile._open_w()
    treantfile.logger.info(
            "task {}, fd {}".format(string, treantfile.handle.fileno()))
    treantfile._close()

    treantfile.add_tags(*["{}_{}".format(string, i) for i in range(1000)])


def test_death_by_1000_pokes(treantfile):
    pool = mp.Pool(processes=4)
    for i in range(12):
        pool.apply_async(pokefile, args=(treantfile.filename,
                                         "run_{}".format(i)))
    pool.close()
    pool.join()

    assert len(treantfile.get_tags()) == 12000
