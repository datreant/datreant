"""
datreant.cli
"""
from __future__ import print_function
import json

import click
import datreant as dtr


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


class AliasedGroup(click.Group):

    def list_commands(self, ctx):
        subcommands = ['init',
                       'show',
                       'draw',
                       'discover',
                       'get',
                       'tags',
                       'categories',
                       'add',
                       'set',
                       'del',
                       'clear',
                       'help']
        return subcommands


@click.command(cls=AliasedGroup)
@click.version_option()
def cli():
    """CLI interface for datreant.

    All subcommands take <dirs> as input, by default using the current
    directory '.' if none given. If '-' given for <dirs>, the subcommand
    will read these from STDIN.

    """
    pass


@cli.command()
@click.pass_context
def help(ctx):
    """Show the CLI help contents"""
    click.echo(ctx.parent.get_help())


def _handle_stdin(dirs):
    if (len(dirs) == 1) and (dirs[0] == '-'):
        try:
            with click.get_text_stream('stdin') as stdin:
                stdin_text = stdin.read()
            dirs = stdin_text.strip().split()
        except:
            pass
    elif len(dirs) == 0:
        dirs = ('.',)

    return dirs


def _set(treant, tags, categories):
    if tags is not None:
        treant.tags = set(tags)
    if categories is not None:
        treant.categories = _parse_categories(categories)


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--tags", '-t', cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def init(dirs, tags=None, categories=None):
    """Make treants from dirs, optionally setting tags/categories"""
    dirs = _handle_stdin(dirs)

    for i in dirs:
        dtr.Treant(i)

    dirs = dtr.Bundle(dirs)
    dirs.map(_set, tags=tags, categories=categories)


def print_treant(treant, detail=False):
    """Print treant with more human readable information"""
    res = dict()
    res['relpath'] = treant.relpath
    res['abspath'] = treant.abspath
    if detail:
        if treant.tags:
            res['tags'] = list(treant.tags)
        if treant.categories:
            res['categories'] = dict(treant.categories)
    return res


@cli.command()
@click.argument("dirs", nargs=-1)
def show(dirs):
    """Show treant metadata contents"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    res = dirs.map(print_treant, detail=True)

    if res:
        click.echo(json.dumps(res, indent=4))


@cli.command()
@click.argument("dirs", nargs=-1)
def draw(dirs):
    """Draw treant filesystem structure"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    for treant in dirs:
        treant.draw()


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--relpath/--abspath", "-r/-a", 'relpath', default=True,
              help="return relative or absolute paths")
@click.option("--tags", '-t', cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
@click.option("--json", '-j', 'json_', is_flag=True,
              help="return results in JSON format")
def discover(dirs, relpath, tags, categories, json_):
    """Search folder for treants and list them"""
    dirs = _handle_stdin(dirs)

    dirs = dtr.Bundle([dtr.discover(i) for i in dirs])

    res = _get(dirs, tags, categories)

    if json_:
        res = res.map(print_treant, detail=True)
        click.echo(json.dumps(res, indent=4))
    else:
        paths = res.relpaths if relpath else res.abspaths
        click.echo("\n".join(paths))


def _get(bundle, tags, categories):
    if tags is None:
        tags = []

    if categories is None:
        categories = {}
    else:
        categories = _parse_categories(categories)

    return bundle.get(*tags, **categories)


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--relpath/--abspath", 'relpath', default=True,
              help="return relative or absolute paths")
@click.option("--json", '-j', 'json_', is_flag=True,
              help="return results in JSON format")
@click.option("--tags", '-t', cls=OptionEatAll,
              help="list of tags or search string")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def get(dirs, relpath, json_, tags, categories):
    """Filter treants that strictly match tags and/or categories

    """
    dirs = _handle_stdin(dirs)

    dirs = dtr.Bundle(dirs, ignore=True)

    res = _get(dirs, tags, categories)

    if json_:
        res = res.map(print_treant, detail=True)
        if res:
            click.echo(json.dumps(res, indent=4))
    else:
        click.echo("\n".join(res.abspaths))


@cli.command()
@click.argument("dirs", metavar='<dir(s)>', nargs=-1)
@click.option("--all/--any", 'all_', default=True,
              help="results for tags present in ALL or ANY treants")
@click.option("--json", '-j', 'json_', is_flag=True,
              help="return results in JSON format")
@click.option("--has", cls=OptionEatAll,
              help=("list of tags; for each returns true if present,"
                    " false if absent"))
def tags(dirs, all_, json_, has):
    """Show treant tags"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)

    if all_:
        res = list(dirs.tags.all)
    else:
        res = list(dirs.tags.any)

    if has:
        res = [tag in res for tag in has]

    if res:
        if json_:
            click.echo(json.dumps(res, indent=4))
        else:
            for i in res:
                click.echo(i)


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--all/--any", 'all_', default=True,
              help="results for keys present in ANY treants")
@click.option("--json", '-j', 'json_', is_flag=True,
              help="results in JSON format")
@click.option("--get", cls=OptionEatAll, help="get values for keys")
def categories(dirs, all_, json_, get):
    """Show treant categories"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)

    if get:
        res = dirs.categories[set(get)]
        keys = get
    else:
        if all_:
            res = dict(dirs.categories.all)
        else:
            res = dict(dirs.categories.any)
        keys = res.keys()

    if res:
        if json_:
            click.echo(json.dumps(res, indent=4))
        else:
            for i, key in enumerate(keys):
                for value in res[key]:
                    click.echo("{} : {}".format(key, value))

                if i < len(res) - 1:
                    click.echo("")


def _parse_categories(catstrings):
    return {key: value for key, value in
            (c.split(":") for c in catstrings)}


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--tags", '-t', cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def add(dirs, tags, categories):
    """Add or update treant tags and/or categories"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    if tags is not None:
        dirs.tags.add(tags)
    if categories is not None:
        dirs.categories.add(_parse_categories(categories))


@cli.command(name='set')
@click.argument("dirs", nargs=-1)
@click.option("--tags", '-t', cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a colon ':'",
)
def set_(dirs, tags, categories):
    """Set/replace all treant tags and/or categories"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    dirs.map(_set, tags=tags, categories=categories)


@cli.command(name='del')
@click.argument("dirs", nargs=-1)
@click.option("--tags", '-t', cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    '-c',
    cls=OptionEatAll,
    help="list of category keys"
)
def del_(dirs, tags, categories):
    """Delete treant tags and/or categories"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    for treant in dirs:
        if tags is not None:
            treant.tags.remove(*tags)
        if categories is not None:
            treant.categories.remove(*categories)


@cli.command()
@click.argument("dirs", nargs=-1)
@click.option("--tags", '-t', is_flag=True, help="remove all tags")
@click.option("--categories", '-c', is_flag=True, help="remove all categories")
def clear(dirs, tags, categories):
    """Clear treant tags and/or categories"""
    dirs = _handle_stdin(dirs)

    # get all existing Treants among dirs
    dirs = dtr.Bundle(dirs, ignore=True)
    for treant in dirs:
        if tags:
            treant.tags.clear()
        if categories:
            treant.categories.clear()
