from pyparsing import (CaselessLiteral, Word, quotedString,
                       removeQuotes, operatorPrecedence, opAssoc, stringEnd,
                       ParseException)

__all__ = ['parse_selection']


class UnaryOperation(object):
    def __init__(self, t):
        self.op, self.a = t[0]


class BinaryOperation(object):
    def __init__(self, t):
        self.op = t[0][1]
        self.operands = t[0][0::2]


class SearchAnd(BinaryOperation):
    def generate_tag_expr(self):
        return list(oper.generate_tag_expr() for oper in self.operands)

    def __repr__(self):
        return "AND:(%s)" % (",".join(str(oper) for oper in self.operands))


class SearchOr(BinaryOperation):
    def generate_tag_expr(self):
        return tuple(oper.generate_tag_expr() for oper in self.operands)

    def __repr__(self):
        return "OR:(%s)" % (",".join(str(oper) for oper in self.operands))


class SearchNot(UnaryOperation):
    def generate_tag_expr(self):
        exps = self.a.generate_tag_expr()
        if isinstance(exps, list):
            return [{e} for e in exps]
        elif isinstance(exps, tuple):
            return tuple({e} for e in exps)
        else:
            return {exps}

    def __repr__(self):
        return "NOT:(%s)" % str(self.a)


class SearchTerm(object):
    def __init__(self, tokens):
        self.term = tokens[0]

    def generate_tag_expr(self):
        return "{}".format(self.term)

    def __repr__(self):
        return self.term


# define the grammar
and_ = CaselessLiteral("and")
or_ = CaselessLiteral("or")
not_ = CaselessLiteral("not")
allowed_chars = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'*+,-./:;<=>?@[]^_`{|}~'
# first remove matching strings and then parse for printable characters
searchTerm = quotedString.setParseAction(removeQuotes) | Word(allowed_chars)
searchTerm.setParseAction(SearchTerm)
searchExpr = operatorPrecedence(searchTerm, [
    (not_, 1, opAssoc.RIGHT, SearchNot),
    (and_, 2, opAssoc.LEFT, SearchAnd),
    (or_, 2, opAssoc.LEFT, SearchOr),
])
Parser = searchExpr + stringEnd


def parse_selection(sel):
    """
    Parse a tag selection string and convert it to the default list/tuple/set
    tag selections. If the selection can't be parsed the original selection
    is returned.

    Parameters
    ----------
    sel : str
        selection string

    Returns
    -------
    list/tuple/set representation used to filter tags

    Example
    -------
    >>> parse_selection('food and drinks')
    ['food', 'drinks']
    >>> parse_selection('free beer')
    'free beer'
    """
    try:
        return Parser.parseString(sel)[0].generate_tag_expr()
    except ParseException:
        return sel
