#! /usr/bin/env python3

# Run usearch -fastq_filter to remove low quality reads

# John Conery
# University of Oregon
# 2014-11-25

from FASTQ import *

import sqlite3
import argparse
import os
import os.path
import sys

###
# Initialize the directory where intermediate work products will be stored.

def init_workspace(args):
	dirname = args.workspace
	
	if not os.path.isdir(dirname):
		os.mkdir(dirname)
	else:
		for root, dirs, files in os.walk(dirname, topdown=False):
			for name in files:
				os.remove(os.path.join(root, name))
			for name in dirs:
				os.rmdir(os.path.join(root, name))

###
# Fetch sample names and the corresponding file names
	
fetch_samples = 'SELECT sample_id, name, r1_file, r2_file FROM samples'

def fetch_sample_list(db, args):
	"Return a list of sample names and the corresponding FASTQ file names."
	sql = fetch_samples
	if isinstance(args.sample, str):
		sql += ' WHERE name = "{}"'.format(args.sample)
	return db.execute(sql).fetchall()

###
# Run usearch on two FASTQ files, merge the results into two new files

def filter_pair(args, sid, fn1, fn2):
	tmp1 = os.path.join(args.workspace, "R1.fastq")
	tmp2 = os.path.join(args.workspace, "R2.fastq")
	
	if not (run_usearch(fn1, tmp1, args) and run_usearch(fn2, tmp2, args)):
		print('usearch failed', fn1, fn2)
		return
		
	f1 = open(tmp1)
	d = dict()
	x = FASTQ.read(f1)
	while x is not None:
		d[x.unique_defline()] = x
		x = FASTQ.read(f1)
	f1.close()
	
	of1 = open(os.path.join(args.workspace, fn1), "w")
	of2 = open(os.path.join(args.workspace, fn2), "w")
	
	f2 = open(tmp2)
	x = FASTQ.read(f2)
	while x is not None:
		if x.unique_defline() in d:
			print(d[x.unique_defline()], file=of1)
			print(x, file=of2)
		x = FASTQ.read(f2)
	f2.close()
	of1.close()
	of2.close()

	
def run_usearch(infile, outfile, args):
	infile = os.path.join(args.directory, infile)
	cmnd = 'usearch -fastq_filter {input} -fastqout {output} -fastq_maxee {ee}'.format(input=infile, output=outfile, ee=args.ee)
	print(cmnd)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
	db.commit()
	return os.system(cmnd) == 0
	

###
# Import the assembled sequences

def import_results(db, args, sid):
	pass
	# file = open(os.path.join(args.workspace, merge_file_pattern.format(sid)))
	# defline = file.readline()
	# count = 0
	# limit = int(args.limit) if args.limit else None
	# while len(defline) > 0:
	# 	# defline = ':'.join(defline[1:].split(':')[:-1])
	# 	defline = FASTQ.parse_defline(defline)
	# 	sequence = file.readline()
	# 	db.execute(insert_sequence, (sid, defline, sequence.strip()))
	# 	count += 1
	# 	if args.limit and count >= limit:
	# 		break
	# 	defline = file.readline()

###
# Top level function: initialize the workspace directory, get sample parameters from
# the database, run usearch for all specified samples

def filter_fastq_files(db, args):
	init_workspace(args)
	for sid, sname, fn1, fn2 in fetch_sample_list(db, args):
		filter_pair(args, sid, fn1, fn2)

###
# Set up command line arguments

def init_api():
	parser = argparse.ArgumentParser(
		description="""Run usearch to filter low quality sequences, write new FASTQ files."""
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-w', '--workspace', help='working directory', default='filtered')
	parser.add_argument('-d', '--directory', help='directory with input FASTQ files', default='data')
	parser.add_argument('-e', '--ee', help='max expected error cutoff', type=float, default=0.3)
	parser.add_argument('-s', '--sample', metavar='id', required=False, help='sample to filter')
	parser.add_argument('-a', '--all', action='store_true', help='filter all samples')
	return parser.parse_args()

###
# Parse the command line arguments, call the top level function...

if __name__ == "__main__":
	args = init_api()
	
	if not (args.sample or args.all):
		argparse.ArgumentParser.exit(1, 'Specify a sample name or --all for all samples')
			
	db = sqlite3.connect(args.dbname)
		
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	filter_fastq_files(db, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )
	
	db.commit()
