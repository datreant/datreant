"""Tests for syncing directories using rsync

"""

import datreant.core as dtr
import pytest
import os
import py.path

def test_sync_locally(tmpdir):
    '''Test that syncronization works locally'''
    with tmpdir.as_cwd():
        sequoia = dtr.Tree("sequoia").makedirs()
        sequoia2 = dtr.Tree("sequoia2").makedirs()
        py.path.local('sequoia/hello.txt').write('hello')
        
        assert os.path.exists('sequoia/hello.txt')
        
        sequoia.sync(sequoia2)
        assert os.path.exists('sequoia2/hello.txt')
        

def test_rsync_tree(tmpdir):

    # We want to sync this tree to another directory
    with tmpdir.as_cwd():
        # Let's assume we have an existing tree
        dtr.Tree("sequoia").makedirs()

        # We create a new tree
        t = dtr.Tree('sequoia2')

        # Options
        t.ssync('sequoia',
               mode='download',
               compress=True,
               backup=True,
               dry=True,
               include=["*.txt"],
               exclude=[])

        # Upload content to sequoia
        t.sync('sequoia', mode='upload')

        # Download content from sequoia
        t.sync('sequoia', mode='download')

        # Remote sync
        t.sync('sequoia2', host='127.0.0.1', user='xxx', password='xxx')
        t.sync('sequoia2', host='127.0.0.1', user='xxx', key='./cloud.pem')
