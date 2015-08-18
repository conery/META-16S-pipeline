#! /usr/bin/env python3

# Detailed comparison of how individual sequences were allocated to clusters
# by different 16S rRNA pipelines

# John Conery
# University of Oregon
# 2014-12-13

# Usage:
#
#	compare_members.py CF1 CF2
#

import argparse
import sqlite3
import os.path
import sys
import re

from FASTA import *

# Global variables...

seqmap1 = { }		# seqmap[d] is the DNA sequence associated with defline d
seqmap2 = { }

###
# UCLUST adds new fields to a defline to show how a sequence was handled.  This
# function extracts the original sequence ID and the info about what happened to
# a sequence (create new OTU, add to existing OTU, chimera)

def parse_defline(s):
	x = s.defline().split(';')		# UCLUST fields separated by semicolon
	return x[0].strip('>'), re.match(r'up=(\w+)',x[-1]).group(1)

###
# Top level...

def compare_members(args):
	cf1 = FASTAReader(args.cf1)
	cf2 = FASTAReader(args.cf2)
	seq1 = cf1.read()
	seq2 = cf2.read()
	while seq1 is not None and seq2 is not None:
		def1, fate1 = parse_defline(seq1)
		def2, fate2 = parse_defline(seq2)
		if def1 != def2:
			raise Exception('Files out of sync: comparing {} with {}'.format(def1,def2))
		seqmap1[def1] = seq1.sequence
		seqmap2[def2] = seq2.sequence
		seq1 = cf1.read()
		seq2 = cf2.read()

###
# Parse command line arguments...

def init_api():
	parser = argparse.ArgumentParser(
		description="...",
		epilog="..."
	)
	parser.add_argument('cf1', help='FASTA file with clusters from baseline pipeline')
	parser.add_argument('cf2', help='FASTA file with clusters from best hit pipeline')
	return parser.parse_args()

###
# Top level....
	
if __name__ == "__main__":
	args = init_api()
	compare_members(args)