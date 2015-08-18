#! /usr/bin/env python3

# Summarize results in a PIP/NGS database.

# John Conery
# University of Oregon
# 2014-00-24

import sqlite3
import argparse
from datetime import *

###
# Read a set of log messages from a file
# TO DO: file should allow either bar or tab as separator

def read_times_from_file(fn):
	with open(fn) as f:
		recs = list(f.readlines())
		starts = filter(lambda x: x.find('.py') > 0, recs)
	return list(map(lambda x: x.split('|'), starts))

###
# Top level:  print execution times of each stage
# TO DO: add start/end entries to log

def print_execution_times(args):
	if args.filename:
		recs = read_times_from_file(args.filename)
	else:
		recs = fetch_times_from_db(args.dbname)
	start = datetime.strptime(recs[0][0], '%Y-%m-%d %H:%M:%S')
	script = recs[0][1]
	for x in recs[1:]:
		end = datetime.strptime(x[0], '%Y-%m-%d %H:%M:%S')
		dt = end - start
		print(dt, script)
		start = end
		script = x[1]
	
###
# Top level:  print stats about data tables

tables = ['reads', 'panda', 'uniq', 'seeds', 'chimeras', 'otu', 'taxonomy']
bcids = ['reads', 'panda', 'uniq', 'otu']

def print_counts(args):
	db = sqlite3.connect(args.dbname)
	if args.ungrouped:
		print_count_matrix(db)
	else:
		print_table_sizes(db)
		
def print_table_sizes(db):
	for t in tables:
		n = db.execute('SELECT count(*) FROM {}'.format(t)).fetchall()[0][0]
		print('{:10s} {:8d}'.format(t,n))
	
### 
# Set up command line arguments

def init_api():
	parser = argparse.ArgumentParser(
		description="Summarize results in a PIP/NGS database",
	)
	parser.add_argument('-d', '--dbname', help='use data from a SQLite database')
	parser.add_argument('-f', '--filename', help='use data from a text file')
	parser.add_argument('-t', '--times', action='store_true', help='print execution times of each stage')
	parser.add_argument('-c', '--counts', action='store_true', help='print the number of items in data tables')
	parser.add_argument('-u', '--ungrouped', action='store_true', help='print stats separately for each experiment')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	
	if not (args.dbname or args.filename):
		print('specify --dbname or --filename')
		exit(1)
	
	if not (args.times or args.counts):
		print("specify --times or --counts")
		exit(1)

	if args.times:
		print_execution_times(args)
	elif args.counts:
		print_counts(args)
