#! /usr/bin/env python3

# Run usearch to filter chimeras using a reference database

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

from common import *

import sqlite3
import argparse
import os
import os.path
import re
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'chimeras'
result_file = 'results.uc'
input_file = 'clusters.fasta'

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
# Initialize (or re-initialize) the table that will hold the output from usearch

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)
create_table = 'CREATE TABLE {tbl} ( cluster_id INTEGER NOT NULL REFERENCES panda, chimeric CHAR(1) )'.format(tbl = table_name)
	
def prepare_table(db, args):
	"""
	Return True if the results table is initialized and ready to accept values.  If the
	table exists already don't overwrite it unless --force was specified on the command line.
	"""
	if db.execute(find_table).fetchall():
		if args.force:
			db.execute(drop_table)
		else:
			return False
	db.execute(create_table)	
	return True	 

###
# Print the sequences to check; they're the centroids in the clusters saved in the clusters table

fetch_clusters = 'SELECT name, sequence FROM clusters'

def print_sequences(db, args):
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', fetch_clusters))
	ff = open(os.path.join(args.workspace, input_file), 'w')
	for name, sequence in db.execute(fetch_clusters):
		print('>{}'.format(name), file=ff)
		print(sequence, file=ff)	
	ff.close()

###
# Run the app

def run_uchime_ref(args):
	cmnd = 'usearch -uchime_ref '
	cmnd += os.path.join(args.workspace, input_file)
	cmnd += ' -db ' + args.reference
	cmnd += ' -strand plus'
	cmnd += ' -uchimeout ' + os.path.join(args.workspace, result_file)
	print(cmnd)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
	db.commit()
	res = os.system(cmnd)
	
###
# Parse the output file, save references to chimeric sequences in a new table.

insert_record = 'INSERT INTO {} (cluster_id, chimeric) VALUES (?,?)'.format(table_name)

def import_results(db, args):
	cmap = make_cmap(db)
	for line in open(os.path.join(args.workspace, result_file)):
		res = line.split('\t')
		defspec = res[1]
		chimeric = res[-1].strip()
		# m = re.search(r'^(.*);', defspec)
		# defline = FASTQ.parse_defline(m.group(1))
		# db.execute(insert_record, (dm[defline], chimeric))
		# db.execute(insert_record, (cmap[m.group(1)], chimeric))
		db.execute(insert_record, (cmap[defspec], chimeric))

def make_cmap(db):
	cmap = { }
	for cid, name in db.execute('SELECT cluster_id, name FROM clusters'):
		cmap[name] = cid
	return cmap

###
# Top level function: initialize the workspace directory, run the app

def filter_chimeras(db, args):
	init_workspace(args)
	print_sequences(db, args)
	run_uchime_ref(args)
	import_results(db, args)

###
# Parse the command line arguments, call the top level function...

def init_api():
	parser = argparse.ArgumentParser(
		description="Run usearch to filter chimeras using a reference sequence",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
	parser.add_argument('-w', '--workspace', help='working directory', default='chimeras')
	parser.add_argument('-d', '--directory', help='directory with input FASTQ files', default='seeds')
	parser.add_argument('-r', '--reference', help='FASTA file containing reference sequences', default=sys.path[0]+'/resources/gold.fa')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	
	db = sqlite3.connect(args.dbname)
		
	if not prepare_table(db, args):
		argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')

	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	filter_chimeras(db, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

	db.commit()
