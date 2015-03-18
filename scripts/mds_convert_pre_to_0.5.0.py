#!/usr/bin/env python2

"""This script will convert state files from pre-release format to 0.5.0 format.
   
The only change is in state file naming: the uuid of the file is included in
the filename.

"""

import mdsynthesis as mds
import os
import tables
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert state file from pre-release format to 0.5.0 format.")
    parser.add_argument("statefile", metavar="STATEFILE", nargs='+', help='statefile(s) to convert')
    args = parser.parse_args()

    for statefile in args.statefile:

        f = tables.open_file(statefile, 'r')
        table = f.get_node('/', 'meta')

        containertype = table.cols.containertype[0]
        uuid = table.cols.uuid[0]
        f.close()

        newfile = os.path.join(os.path.dirname(statefile), '{}.{}.h5'.format(containertype, uuid))
        os.rename(statefile, newfile)

