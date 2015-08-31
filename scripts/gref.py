#! /usr/bin/env python3

# Search a FASTA or FASTQ file for a sequence that matches a pattern.  The
# pattern can be in the defline or the sequence letters.  Inspired by the 
# 'gref' script in the SEALS package.

# John Conery
# University of Oregon
# 2015-08-22

# TBD: 2nd arg optional; if not supplied, read from stdin

import argparse
import os
import sys
import re

from FASTA import *
from FASTQ import *

# When searching a DNA file the script can use a regular expressions based on IUPAC
# ambiguity letter to match the sequence

iupac_map = {
	'A' : r'A',	 
	'C' : r'C',	 
	'G' : r'G',	 
	'T' : r'T',
	'R' : r'[AG]', 
	'Y' : r'[CT]',
	'S' : r'[GC]',
	'W' : r'[AT]',
	'K' : r'[GT]',
	'M' : r'[AC]',
	'B' : r'[CGT]',
	'D' : r'[AGT]',
	'H' : r'[ACT]',
	'V' : r'[ACG]',
	'N' : r'[ACGT]'
}

def iupac_regexp(pat):
    'Create a regular expression for a sequence, replacing IUPAC letters with groups'
    return ''.join(map(lambda ch: iupac_map[ch], pat))
    
complement = {
    'A' : 'T', 
    'T' : 'A', 
    'C' : 'G',
    'G' : 'C',
	'R' : r'[TC]', 
	'Y' : r'[GA]',
	'S' : r'[AT]',
	'W' : r'[GC]',
	'K' : r'[CA]',
	'M' : r'[TG]',
	'B' : r'.',
	'D' : r'.',
	'H' : r'.',
	'V' : r'.',
	'N' : r'.'
}

def reverse_complement(pat):
    return ''.join(map(lambda ch: complement[ch], reversed(pat)))

###
# Parse command line arguments...

def init_api():
    parser = argparse.ArgumentParser(
        description="""Find sequences that match a pattern in FASTA or FASTQ files.  If the --seq option
is used the program looks for the pattern in sequence characters, otherwise it looks for the pattern in
deflines.  Patterns are specified using PERL regular expression syntax.         
If the inputs are named files the program will infer the sequence type from the filename extension 
(.fa, .fn, or .fasta for FASTA files, .fastq for FASTQ files).  If the input is coming from stdin the 
sequence type must be defined with the --fasta or --fastq option.  
"""
)
    parser.add_argument('pattern', help='pattern to find')
    parser.add_argument('files', nargs='*', help='name(s) of sequence file(s)')
    parser.add_argument('-s', '--sequence', action='store_true', help='match sequence characters')
    parser.add_argument('-i', '--iupac', action='store_true', help='pattern contains IUPAC ambiguity letters')
    parser.add_argument('-v', action='store_true', help='print sequences that do not match the pattern')
    parser.add_argument('-r', '--reverse', action='store_true', help='use the reverse complement of the pattern')
    parser.add_argument('--fasta', action='store_true', help='input stream is in FASTA format')
    parser.add_argument('--fastq', action='store_true', help='input stream is in FASTQ format')
    return parser.parse_args()

reader_type = {
    '.fa'    : FASTAReader,
    '.fn'    : FASTAReader,
    '.fasta' : FASTAReader,
    '.fastq' : FASTQReader,
}

# The string %s inserted into this pattern will be printed in red (or the terminal's ANSI color #1)

ESC = chr(27)

highlight_region = ESC + '[1;31m' + '%s' + ESC + '[0m'

###
# Check the command line arguments

def validate(args):
    # if there are no file names we need either --fasta or --fastq
    if len(args.files) == 0:
        if bool(args.fasta) == bool(args.fastq):
            argparse.ArgumentParser.exit(1, 'specify either --fasta or --fastq when reading from stdin')
    # all file names need to end in a recognized extension
    for fn in args.files:
        if os.path.splitext(fn)[1] not in reader_type:
            argparse.ArgumentParser.exit(1, 'invalid filename extension: ' + fn)

###
# Process a single file

def scan_file(file, pattern, tty, args):
    for seq in file:
        src = seq.sequence() if args.sequence else seq.defline()
        match = re.search(pattern, src)
        if match is not None and not args.v:
            res = repr(seq)
            hit = match.group()
            if print_to_tty:
                res = re.sub(hit, highlight_region % hit, res) 
            print(res)
        elif match is None and args.v:
            print(repr(seq))
      
###
# Main program:
    
if __name__ == "__main__":
    args = init_api()
    validate(args)
    
    if args.iupac:
        pattern = iupac_regexp(args.pattern) 
    elif args.reverse:
        pattern = reverse_complement(args.pattern)
    else:
        pattern = args.pattern
        
    print_to_tty = os.isatty(sys.stdout.fileno())
    
    if len(args.files) == 0:
        reader = FASTAReader(sys.stdin) if args.fasta else FASTQReader(sys.stdin)
        scan_file(reader, pattern, print_to_tty, args)
    else:
        for fn in args.files:
            ext = os.path.splitext(fn)[1]
            cls = reader_type.get(ext)
            scan_file(cls(fn), pattern, print_to_tty, args)
    


