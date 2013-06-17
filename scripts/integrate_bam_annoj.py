#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
import string
import sys
import argparse

from pypeline import Pypeline, ProcessPypeStep as PPS, PythonPypeStep as PyPS, parse_steplist

usage = "%prog aligner [options]"
description = """Takes a BAM alignment file, insert into a AnnoJ database 
and creates a fetcher.
"""

epilog = """Note: it is advised to leave the --*-args arguments unchanged
unless you really know what you're doing."""

parser = argparse.ArgumentParser(description=description,epilog=epilog)
subparsers = parser.add_subparsers(dest='aligner_name')

parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--dbhost',dest='dbhost',action='store',help='host of AnnoJ database')
parent_parser.add_argument('--dbname',dest='dbname',action='store',help='name of AnnoJ database')
parent_parser.add_argument('--dbuser',dest='dbuser',action='store',help='user name of AnnoJ database')
parent_parser.add_argument('--dbpass',dest='dbpass',action='store',help='password of AnnoJ database')
parent_parser.add_argument('--dbtable',dest='dbtable',action='store',required=True,help='table name to use in the AnnoJ database')
parent_parser.add_argument('--bamfile',dest='bamfile',action='store',required=True,help='name of alignment file in BAM format')
parent_parser.add_argument('--fetcher',dest='fetcher',action='store',required=True,help='file name of fetcher')
parent_parser.add_argument('--title',dest='title',action='store',help='title of AnnoJ browser track')
parent_parser.add_argument('--info',dest='info',action='store',help='info line of AnnoJ browser track')

parser_bowtie = subparsers.add_parser('bowtie',parents=[parent_parser])
parser_bowtie2 = subparsers.add_parser('bowtie2',parents=[parent_parser])
parser_tophat = subparsers.add_parser('tophat',parents=[parent_parser])
parser_tophat.add_argument('--pe',dest='pe',action='store_true',help='whether the BAM file represents paired-end runs')
parser_tophat2 = subparsers.add_parser('tophat2',parents=[parent_parser])
parser_tophat2.add_argument('--pe',dest='pe',action='store_true',help='whether the BAM file represents paired-end runs')

if __name__ == '__main__' :

    # parse command line arguments
    args = parser.parse_args(sys.argv[1:])
    if (args.title==None):
        args.title = args.dbtable
    if (args.info==None):
        args.info = args.dbtable


