"""
datreant.cli
"""
from __future__ import print_function

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
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


@click.group()
@click.version_option()
def cli():
    """CLI interface for datreant filesystem databases"""
    pass


@cli.command()
@click.argument("folder", default=".")
def init(folder):
    """turn folder into a treant"""
    dtr.Treant(folder)


def print_treant(treant, verbose=False):
    """print treant with more human readable information"""
    print("abspath: ", treant)
    if verbose:
        print("tags: ", treant.tags)
        print("categories: ", treant.categories)


@cli.command()
@click.argument("folder", default=".")
def show(folder):
    """show content of treant"""
    tree = dtr.Tree(folder)
    if tree[".datreant"].exists:
        treant = dtr.Treant(folder)
        print_treant(treant, verbose=True)
    else:
        print("not a treant")


@cli.command()
@click.option("--tags", cls=OptionEatAll, help="list of tags")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a double point ':'",
)
@click.argument("folder", default=".")
def update(tags, categories, folder):
    """update tags and categories of treant"""
    treant = dtr.Treant(folder)
    if tags is not None:
        treant.tags = set(tags)
    if categories is not None:
        for key, value in (c.split(":") for c in categories):
            treant.categories[key] = value


@cli.command()
@click.argument("folder", default=".")
@click.option("--tags", cls=OptionEatAll, help="list of tags or search string")
@click.option(
    "--categories",
    cls=OptionEatAll,
    help="list of categories as key-value pairs separated by a double point ':'",
)
@click.option("--verbose", is_flag=True, help="list tags and categories as well")
def search(tags, folder, categories, verbose):
    """search folder for treants and list them"""
    bundle = dtr.discover(folder)

    if tags is not None:
        if len(tags) != 1:
            tags = set(tags)
        bundle = bundle[bundle.tags[tags]]

    if categories is not None:
        categories = {k: v for k, v in (c.split(":") for c in categories)}
        groupby = bundle.categories.groupby(list(categories.keys()))
        bundle = groupby[set(categories.values())]

    for treant in bundle:
        if verbose:
            print_treant(treant, verbose=True)
        else:
            print(treant.abspath)


# for easier debugging
if __name__ == "__main__":
    cli()
