#! /usr/bin/env python3


# Print stats about sequences in a file

# John Conery
# University of Oregon
# 2014-12-13

# Usage:
#
#    fstats.py file
#

import argparse

from FASTA import *
from FASTQ import *

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
    
    if args.fn.endswith('.fasta'):
        reader = FASTAReader(args.fn)
    elif args.fn.endswith('.fastq'):
        reader = FASTQReader(args.fn)
    else:
        print('file name must end with .fasta or .fastq')
        exit(1)
        
    stats = { 'n' : 0, 'min': 10000, 'max': 0, 'sum' : 0}
    
    for seq in reader:
        n = len(seq)
        stats['n'] += 1
        stats['sum'] += n
        stats['min'] = min(n, stats['min'])
        stats['max'] = max(n, stats['max'])
    
    print(stats['n'], 'sequences')
    print('shortest: {:7d}'.format(stats['min']))
    print('longest:  {:7d}'.format(stats['max']))
    print('mean:     {:7d}'.format(stats['sum']//stats['n']))

