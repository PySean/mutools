#!/usr/bin/env python
"""
    "reindex.py", by Sean Soderman
    Reindexes an entire directory of BAM files.
"""
import os
import sys
from subprocess import check_output

if len(sys.argv) < 2:
    sys.stderr.write('Usage: {} <bamdir>'.format(sys.argv[0]))
    sys.exit(1)

directory = sys.argv[1]
for dirpath, dirnames, filenames in os.walk(directory):
    for f in filenames:
        if f.endswith('.bam'):
            filepath = os.path.join(dirpath, f)
            try:
                check_output(['samtools', 'index', filepath])
            except subprocess.CalledProcessError as cpe:
                sys.stderr.write("Error: Samtools messed up: {}".format(cpe))
                sys.exit(1)
