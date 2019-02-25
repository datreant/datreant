"""Tests for the chicken licking interface

"""

import subprocess
import pathlib
import pytest

import datreant as dtr


@pytest.fixture()
def in_tmpdir(tmpdir):
    with tmpdir.as_cwd():
        yield


@pytest.mark.parametrize('tags', [False, True])
@pytest.mark.parametrize('categories', [False, True])
def test_init(in_tmpdir, tags, categories):
    cmd = 'dtr init thisdir'
    if tags:
        cmd += ' -t tasty'
    if categories:
        cmd += ' -c flavour:chicken'

    subprocess.run(cmd,
                   shell=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                   check=True)

    loc = pathlib.Path('./thisdir')
    assert loc.exists
    # Check that treant exists before opening
    dtrloc = (loc / '.datreant')
    assert dtrloc.exists

    t = dtr.Treant('thisdir')

    if tags:
        assert 'tasty' in t.tags
    if categories:
        assert t.categories['flavour'] == 'chicken'


@pytest.fixture
def readymade_treant(in_tmpdir):
    t = dtr.Treant('duchamp',
                   tags='art?',
                   categories={'colour': 'white',
                               'material': 'porcelain'},
    )

def test_show(readymade_treant):
    ret = subprocess.run('dtr show duchamp', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True)

    output = ret.stdout.decode()
    assert 'art?' in output
    assert 'colour' in output
    assert 'white' in output
