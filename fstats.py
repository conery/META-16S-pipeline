#! /usr/bin/env python3


# Print stats about sequences in a file

# John Conery
# University of Oregon
# 2014-12-13

# Usage:
#
#    fstats.pt file
#

import argparse

from FASTA import *
from FASTQ import *

def update(stats, seq):
    n = len(seq)
    stats['n'] += 1
    stats['sum'] += n
    stats['min'] = min(n, stats['min'])
    stats['max'] = max(n, stats['max'])

def compute_fasta_stats(fn, stats):
    f = FASTAReader(fn)
    seq = f.read()
    while seq is not None:
        update(stats, seq.sequence())
        seq = f.read()

def compute_fastq_stats(fn, stats):
    f = open(fn)
    seq = FASTQ.read(f)
    while seq is not None:
        update(stats, seq.sequence())
        seq = FASTQ.read(f)

###
# Parse command line arguments...

def init_api():
    parser = argparse.ArgumentParser(
        description="...",
        epilog="..."
    )
    parser.add_argument('fn', help='name of sequence file (.fasta or .fastq)')
    return parser.parse_args()

###
# Top level....
    
if __name__ == "__main__":
    args = init_api()
    
    stats = { 'n' : 0, 'min': 10000, 'max': 0, 'sum' : 0}
    
    if args.fn.endswith('.fasta'):
        compute_fasta_stats(args.fn, stats)
    elif args.fn.endswith('.fastq'):
        compute_fastq_stats(args.fn, stats)
    else:
        print('file name must end with .fasta or .fastq')
        exit(1)
    
    print(stats['n'], 'sequences')
    print('shortest: {:7d}'.format(stats['min']))
    print('longest:  {:7d}'.format(stats['max']))
    print('mean:     {:7d}'.format(stats['sum']//stats['n']))

