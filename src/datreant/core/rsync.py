"""Rsync wrapper"""

import subprocess
import os

def rsync(source, dest, compress=True, backup=False, dry=False, include=None, exclude=None, rsync_path='/usr/bin/rsync'):
    opts = ['-r'] # Recursive
    
    if compress:
        opts.append('-c')
    
    if backup:
        opts.append('-b')
    
    if dry:
        opts.append('--dry-run')
    
    if include:
        opts.extend(['--include ', include])
    
    if exclude:
        opts.extend(['--exclude ', include])
    
    if not os.path.isdir(source):
        raise ValueError("Source should be a directory")
    
    source = os.path.join(source, '') # Add trailing slash

    p = subprocess.Popen([rsync_path] + opts + [source, dest], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise Exception('Syncing error: rsync returned status code {}\nStandard error: {}\nStandard output: {}'.format(p.returncode, err, out))
