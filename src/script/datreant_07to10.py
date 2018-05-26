#! /usr/bin/env python
from glob import glob
from os import path

import json
import datreant


def convert(folder):
    json_file = glob(path.join(folder, 'Treant*'))[0]
    with open(json_file) as fh:
        old = json.load(fh)
    datreant.Treant(folder, categories=old['categories'], tags=old['tags'])


if __name__ == '__main__':
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
