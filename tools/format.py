#! /usr/bin/env python

"""
Code formatting tool.  Runs uncrustify for constraints that get
enforced automatically, then checks for remaining violations.
"""

import argparse
import os
import sys
import glob

SCRIPTDIR = os.path.normpath(os.path.dirname(__file__))
parser = argparse.ArgumentParser(description=__doc__.strip())
parser.add_argument('--check', dest='check', action='store_true',
                    help='Check the code but do not mutate it')


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def check_format(filename):
    infile = open(filename)
    lineno = 1

    def fail(msg):
        print(bcolors.FAIL + msg.format(filename, lineno) + bcolors.ENDC)

    bad = False
    previous = ''
    previous_is_blank = False
    for line in infile:
        line = line.rstrip('\n')
        if len(line) > 80:
            fail('{}:{} Line is over 80 chars.')
            bad = True
        is_blank = len(line) == 0
        if previous_is_blank and line.lstrip(' ') == '}':
            fail('{}:{} Extra newline before ending brace.')
            bad = True
        previous_is_blank = is_blank
        previous = line
        lineno = lineno + 1

    return not bad

if __name__ == '__main__':
    args = parser.parse_args()
    cmd = 'uncrustify -l C -c {}/.uncrustify *.h '.format(SCRIPTDIR)
    if os.path.exists('uncrustify'):
        cmd = './' + cmd
    if args.check:
        cmd += '--check'
    else:
        cmd += '--no-backup'
    if os.system(cmd) and args.check:
        sys.exit(1)
    good = True
    for filename in glob.glob('*.h'):
        good = check_format(filename) and good
    if not good:
        print 'Illegal formatting detected.'
        sys.exit(1)
