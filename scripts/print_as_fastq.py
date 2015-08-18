#! /usr/bin/env python3

# Print the reads in a PIP/NGS database.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

import sqlite3
import argparse

from FASTQ import *

fetch_reads = 'SELECT read_id, read FROM reads'
	
def print_sequences(args):
	db = sqlite3.connect(args.dbname)
	sql = fetch_reads
	if args.limit is not None:
		sql += ' LIMIT {}'.format(args.limit)
	for rid, blob in db.execute(sql).fetchall():
		x = FASTQ(blob)
		x._def = "{}".format(rid)
		x.unpack()
		print(x)
	
### 
# Set up command line arguments

def init_api():
	parser = argparse.ArgumentParser(
		description="Print the sequences in a PIP/NGS database.",
		epilog="If a limit N is specified the first N reads from both ends of paired end read are printed."
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-l', '--limit', required=False, help='number of sequences to print')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	print_sequences(args)
