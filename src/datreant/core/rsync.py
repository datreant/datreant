"""Rsync wrapper"""

import subprocess
import os


def rsync(source, dest, compress=True, backup=False, dry=False, checksum=True,
          include=None, exclude=None, overwrite=False,
          rsync_path='/usr/bin/rsync'):
    """Wrapper function for rsync. There are some minor differences with the
    standard rsync behaviour:

    - The include option will exclude every other path.
    - The exclude statements will take place after the include statements

    Parameters
    ----------
    source : str
        Source directory for the sync
    dest : str
        Dest directory for the sync
    compress : bool
        If True, use gzip compression to reduce the data transfered over
        the network
    backup : bool
        If True, pre-existing files are renamed with a "~" extension before
        they are replaced
    dry: bool
        If True, do a dry-run. Useful for debugging purposes
    checksum: bool
        Perform a checksum to determine if file content has changed.
    include: str or list
        Paths (wildcards are allowed) to be included. If this option is used,
        every other path is excluded
    exclude: str or list
        Paths to be excluded from the copy
    overwrite: bool
        If False, files in `dest` that are newer than files in `source`
        will not be overwritten
    rsync_path: str
        Path where to find the rsync executable
    """
    opts = ['-r']  # Recursive

    if not overwrite:
        opts.append('--ignore-existing')

    if compress:
        opts.append('-z')

    if checksum:
        opts.append('-c')

    if backup:
        opts.append('-b')

    if dry:
        opts.append('--dry-run')

    if include:
        opts.extend(['--include=*/'])

        if isinstance(include, list):
            opts.extend(["--include={}".format(inc) for inc in include])
        elif isinstance(include, str):
            opts.append("--include={}".format(include))

        opts.extend(['--exclude=*'])

    if exclude:
        if isinstance(exclude, list):
            opts.extend(["--exclude={}".format(exc) for exc in exclude])
        elif isinstance(exclude, str):
            opts.append("--exclude={}".format(exclude))

    source = os.path.join(source, '')  # Add trailing slash
    cmd = [rsync_path] + opts + [source, dest]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    if p.returncode != 0:
        raise Exception('\n'.join([
            'Syncing error: rsync returned status code {}',
            'Standard error: {}'
            'Standard output: {}']).format(
            p.returncode, err, out))

    return cmd
