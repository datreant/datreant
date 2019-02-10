"""
datreant.cli
"""
from __future__ import print_function
import json

import click
import datreant as dtr
from datreant.names import TREANTDIR_NAME


class OptionEatAll(click.Option):
    # from
    # https://stackoverflow.com/a/48394004/2207958

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, "nargs, if set, must be -1 not {}".format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = (parser._long_opt.get(name) or
                          parser._short_opt.get(name))
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


@click.group()
@click.version_option()
def cli():
    """CLI interface for treants"""
    pass


@cli.command()
@click.pass_context
def help(ctx):
    """Show the CLI help contents"""
    click.echo(ctx.parent.get_help())


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
def init(dirs):
    """Turn directory into a treant"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    for i in dirs:
        dtr.Treant(i)


def print_treant(treant, verbose=False):
    """Print treant with more human readable information"""
    res = dict()
    res['relpath'] = treant.relpath
    res['abspath'] = treant.abspath
    if verbose:
        if treant.tags:
            res['tags'] = list(treant.tags)
        if treant.categories:
            res['categories'] = dict(treant.categories)
    return res


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
def show(dirs):
    """Show content of treant"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    res = []
    for i in dirs:
        tree = dtr.Tree(i)
        if tree[TREANTDIR_NAME].exists:
            treant = dtr.Treant(i)
            res.append(print_treant(treant, verbose=True))
    click.echo(json.dumps(res, indent=4))


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
def draw(dirs):
    """Draw filesystem structure of treant"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    res = []
    for treant in dirs:
        treant.draw()


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--melt", is_flag=True, help="merge tags from multiple treants")
@click.option("--json", 'json_', is_flag=True, help="print results in JSON format")
def tags(dirs, melt, json_):
    """Show treant tags"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    res = []
    for treant in dirs:
        if treant.tags:
            res.append(list(treant.tags))

    if melt:
        res_set = set()
        for j in res:
            res_set.update(set(j))
        res = list(res_set)

        if res:
            if json_:
                click.echo(json.dumps(res, indent=4))
            else:
                for tag in res:
                    click.echo(tag)
    else:
        if res:
            if json_:
                click.echo(json.dumps(res, indent=4))
            else:
                for tags in res:
                    for tag in tags:
                        click.echo(tag)
                    click.echo("")


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--json", 'json_', is_flag=True, help="print results in JSON format")
@click.option("--get", cls=OptionEatAll, help="get values for given keys")
def categories(dirs, json_, get):
    """Show treant categories"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    res = []
    for treant in dirs:
        if treant.categories:
            if get:
                res.append(dict(treant.categories[set(get)]))
            else:
                res.append(dict(treant.categories))

    if res:
        if json_:
            click.echo(json.dumps(res, indent=4))
        else:
            for cats in res:
                for key, value in cats.items():
                    click.echo("{} : {}".format(key, value))
                click.echo("")


def _parse_categories(catstrings):
    return {key: json.loads(value) for key, value in
            (c.split(":") for c in catstrings)}


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--tags", cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def add(tags, categories, dirs):
    """Update tags and categories of treant"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    for treant in dirs:
        if tags is not None:
            treant.tags.add(tags)
        if categories is not None:
            treant.categories.add(_parse_categories(categories))


@cli.command(name='set')
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--tags", cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def set_(tags, categories, dirs):
    """Replace tags and/or categories as given"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    for treant in dirs:
        if tags is not None:
            treant.tags = set(tags)
        if categories is not None:
            treant.categories = _parse_categories(categories)


@cli.command(name='del')
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--tags", cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of category keys"
)
def del_(tags, categories, dirs):
    """Delete tags and categories of treant"""
    if not dirs:
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs)
    for treant in dirs:
        if tags is not None:
            treant.tags.remove(*tags)
        if categories is not None:
            treant.categories.remove(*categories)


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--tags", cls=OptionEatAll, help="list of tags or search string")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
@click.option("--verbose", '-v', is_flag=True,
              help="list tags and categories as well")
def search(dirs, tags, categories, verbose):
    """Search folder for treants and list them"""
    bundle = dtr.Bundle([dtr.discover(i) for i in dirs])

    if tags is None:
        tags = []

    if categories is None:
        categories = {}
    else:
        categories = _parse_categories(categories)

    res = bundle.get(*tags, **categories)

    if verbose:
        res = print_treant(treant, verbose=True)
    else:
        click.echo("\n".join(res.abspaths))
