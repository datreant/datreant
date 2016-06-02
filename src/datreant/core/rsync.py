"""Rsync wrapper"""

import subprocess
import os


def rsync(source, dest, compress=True, backup=False, dry=False,
          include=None, exclude=None, rsync_path='/usr/bin/rsync'):
    opts = ['-r']  # Recursive

    if compress:
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
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise Exception('\n'.join([
            'Syncing error: rsync returned status code {}',
            'Standard error: {}'
            'Standard output: {}']).format(
            p.returncode, err, out))
