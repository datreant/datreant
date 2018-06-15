#! /usr/bin/env python
import json
import warnings
from glob import glob
from os import path

import datreant


def convert(folder):
    treants = glob(path.join(folder, 'Treant*'))
    if len(treants) == 0:
        warnings.warn("No treant found in folder: {}".format(folder))
        return
    elif len(treants) > 1:
        warnings.warn("Multiple treants found in folder: {}\n"
                      "    Skipping conversion".format(folder))
        return

    json_file = treants[0]
    with open(json_file) as fh:
        old = json.load(fh)
    datreant.Treant(folder, categories=old['categories'], tags=old['tags'])


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Convert pre 1.0 Treant to new Treant",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "directories",
        metavar="DIRECTORY",
        nargs="+",
        help="one or more directories that are being processed")
    args = parser.parse_args()

    for dir in args.directories:
        convert(dir)


if __name__ == '__main__':
    main()
