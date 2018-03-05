"""Parser tests for Treants.

"""

from datreant.selectionparser import parse_selection

import pytest


@pytest.mark.parametrize('expr', ('foo bar', 'foo'))
def test_names(expr):
    assert parse_selection(expr) == expr


@pytest.mark.parametrize("expr, typename", (("foo and bar", list),
                                            ("foo or bar", tuple), ("not foo",
                                                                    set)))
def test_type_operation(expr, typename):
    sel = parse_selection(expr)
    assert isinstance(sel, typename)


@pytest.mark.parametrize('terms', (['foo', 'bar2'], ['a', 'b', 'c'],
                                   ['d', 'e', 'f', 'g']))
def test_and(terms):
    expr = ' and '.join(terms)
    assert parse_selection(expr) == terms


@pytest.mark.parametrize('terms', (('foo', 'bar2'), ('a', 'b', 'c'),
                                   ('d', 'e', 'f', 'g')))
def test_or(terms):
    expr = ' or '.join(terms)
    assert parse_selection(expr) == terms


@pytest.mark.parametrize('terms', (('foo', 'bar2'), ('a', 'b', 'c'),
                                   ('d', 'e', 'f', 'g')))
def test_not(terms):
    expr = 'not ({})'.format(' and '.join(terms))
    assert parse_selection(expr) == [{el} for el in terms]


def test_complex_expr():
    expr = "not food or (leaves and animals)"
    assert parse_selection(expr) == ({'food'}, ['leaves', 'animals'])


@pytest.mark.parametrize(
    'char',
    list(
        '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'*+,-./:;<=>?@[]^_`{|}~'
    ))
def test_possible_characters(char):
    sel = parse_selection("food and a{}".format(char))
    assert isinstance(sel, list)
