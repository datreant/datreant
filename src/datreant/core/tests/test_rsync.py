"""Tests for syncing directories using rsync

"""

import datreant.core as dtr
import pytest
import os
import py.path


def test_sync_local(tmpdir):
    '''Test that syncronization works locally'''
    with tmpdir.as_cwd():
        sequoia = dtr.Tree("sequoia").makedirs()
        sequoia2 = dtr.Tree("sequoia2").makedirs()
        py.path.local('sequoia/hello.txt').write('hello')

        assert os.path.exists('sequoia/hello.txt')

        sequoia.sync(sequoia2)
        assert os.path.exists('sequoia2/hello.txt')


def _create_tree(name, files=[]):
    tree = dtr.Tree(name).makedirs()
    for file_ in files:
        py.path.local(os.path.join(tree.abspath, file_)
                      ).write('hello', ensure=True)

    return tree


def test_overwrite(tmpdir):
    with tmpdir.as_cwd():
        sequoia = _create_tree('sequoia', ['hello.txt', 'data/hello.txt',
                                           'data/world.dat', 'world.dat'])
        sequoia2 = dtr.Tree("sequoia2").makedirs()

        # Upload contents
        sequoia.sync(sequoia2, mode='upload')

        # Change contents
        open(sequoia2['hello.txt'].abspath, 'w').write('newcontent')

        # Upload contents again
        sequoia.sync(sequoia2, mode='upload')

        # Verify that the hello.txt is not overwritten
        assert sequoia2['hello.txt'].read() == 'newcontent'


def test_excludes(tmpdir):
    with tmpdir.as_cwd():
        sequoia = _create_tree('sequoia', ['hello.txt', 'data/hello.txt',
                                           'data/world.dat', 'world.dat'])
        sequoia2 = dtr.Tree("sequoia2").makedirs()
        sequoia3 = dtr.Tree("sequoia3").makedirs()

        sequoia.sync(sequoia2, exclude='*.txt')

        assert os.path.exists('sequoia2/world.dat')
        assert os.path.exists('sequoia2/data/world.dat')

        assert not os.path.exists('sequoia2/hello.txt')
        assert not os.path.exists('sequoia2/data/hello.txt')

        sequoia.sync(sequoia3, exclude=['*.txt', '*.dat'])
        assert not os.path.exists('sequoia3/hello.txt')
        assert not os.path.exists('sequoia3/world.dat')

        assert os.path.exists('sequoia3/data/')


def test_includes(tmpdir):
    with tmpdir.as_cwd():
        sequoia = _create_tree('sequoia', ['hello.txt', 'data/hello.txt',
                                           'data/world.dat', 'world.dat'])
        sequoia2 = dtr.Tree("sequoia2").makedirs()
        sequoia3 = dtr.Tree("sequoia3").makedirs()
        sequoia4 = dtr.Tree("sequoia4").makedirs()

        # Only txt
        sequoia.sync(sequoia2, include='*.txt')

        assert os.path.exists('sequoia2/data/hello.txt')
        assert os.path.exists('sequoia2/hello.txt')

        assert not os.path.exists('sequoia2/world.dat')
        assert not os.path.exists('sequoia2/data/world.dat')

        # Only txt and dat
        sequoia.sync(sequoia3, include=['*.txt', '*.dat'])

        assert os.path.exists('sequoia3/data/hello.txt')
        assert os.path.exists('sequoia3/hello.txt')

        assert os.path.exists('sequoia3/world.dat')
        assert os.path.exists('sequoia3/data/world.dat')

        # We can also test include/excludes at the same time
        sequoia.sync(sequoia4, exclude='*.txt', include=['data/*'])

        assert os.path.exists('sequoia4/data/world.dat')
        assert os.path.exists('sequoia4/data/hello.txt')
        assert not os.path.exists('sequoia4/hello.txt')
        assert not os.path.exists('sequoia4/world.dat')


# def test_sync_remote(tmpdir):
#     remote_path = 'user@host:/tmp/sequoia'
#     with tmpdir.as_cwd():
#         sequoia = dtr.Tree("sequoia").makedirs()
#         py.path.local('sequoia/hello.txt').write('hello')
#
#         sequoia2 = dtr.Tree("sequoia2").makedirs()
#
#         # Upload
#         sequoia.sync(remote_path)
#
#         # Download
#         sequoia2.sync(remote_path, mode='download')
#
#         assert os.path.exists('sequoia2/hello.txt')


# def test_rsync_tree(tmpdir):
#
#     # We want to sync this tree to another directory
#     with tmpdir.as_cwd():
#         # Let's assume we have an existing tree
#         dtr.Tree("sequoia").makedirs()
#
#         # We create a new tree
#         t = dtr.Tree('sequoia2')
#
#         # Options
#         t.sync('sequoia',
#                mode='download',
#                compress=True,
#                backup=True,
#                dry=True,
#                include=["*.txt"],
#                exclude=[])
#
#         # Upload content to sequoia
#         t.sync('sequoia', mode='upload')
#
#         # Download content from sequoia
#         t.sync('sequoia', mode='download')
#
#         # Remote sync
#         t.sync('sequoia2', host='127.0.0.1', user='xxx', password='xxx')
#         t.sync('sequoia2', host='127.0.0.1', user='xxx', key='./cloud.pem')
