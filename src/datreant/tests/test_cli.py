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
    t = dtr.Treant('duchamp/fountain',
                   tags=['art?', 'duchamp'],
                   categories={'colour': 'white',
                               'material': 'porcelain'},
    )

@pytest.fixture
def readymades(readymade_treant):
    t2 = dtr.Treant('duchamp')
    t3 = dtr.Treant('duchamp/wheel',
                    tags=['duchamp'],
                    categories={'material': 'wood'})
    t4 = dtr.Treant('manray')
    t5 = dtr.Treant('manray/gift',
                    tags=['manray'],
                    categories={'material': 'metal',
                                'spiky': True})
    t6 = dtr.Treant('hausmann')
    t7 = dtr.Treant('hausmann/head',
                    tags=['hausmann'],
                    categories={'material': 'wood',
                                'spiky': False})


def test_show(readymade_treant):
    ret = subprocess.run('dtr show duchamp/fountain', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True)

    output = ret.stdout.decode()
    assert 'art?' in output
    assert 'colour' in output
    assert 'white' in output


def test_discover(readymades):
    ret = subprocess.run('dtr discover', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    output = ret.stdout.decode()

    print(output)
    items = [val for val in output.split('\n')
             if val]  # ignore blank lines
    assert len(items) == 7


def test_get(readymades):
    ret = subprocess.run('dtr get -c material:wood', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    # should return 2 results, 'duchamp/wheel' and 'hausmann/head'
    output = ret.stdout.decode()
    print('#' + output + '#')
    assert len(output.split('\n')) == 2
    assert 'wheel' in output
    assert 'head' in output


def test_tags(readymade_treant):
    ret = subprocess.run('dtr tags duchamp/fountain', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    output = ret.stdout.decode()

    assert 'art?' in output
    assert 'duchamp' in output


def test_categories(readymade_treant):
    ret = subprocess.run('dtr categories duchamp/fountain', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    output = ret.stdout.decode()

    assert 'colour : white' in output


def test_add_tag(readymade_treant):
    ret = subprocess.run('dtr add duchamp/fountain -t dada',
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    t = dtr.Treant('duchamp/fountain')

    assert 'dada' in t.tags


def test_add_category(readymade_treant):
    ret = subprocess.run('dtr add duchamp/fountain -c year:1917',
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True,
    )

    t = dtr.Treant('duchamp/fountain')

    assert t.categories['year'] == 1917
