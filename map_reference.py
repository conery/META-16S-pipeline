#! /usr/bin/env python3

# Use usearch_local to find matches for the unique sequences in the reference database.
# Import the results into a new table named 'hits'

# John Conery
# University of Oregon
# 2014-10-30

from common import *

import sqlite3
import argparse
import os
import os.path
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'hits'
input_file = 'uniq.fasta'
output_file = 'hits.tsv'

# ### Deprecated -- table used by version that imported all hits into DB ##
# # The tuples in this last are the field names passed to usearch_local and the corresponding
# # column names for that value in the hits table
#
# hit_fields = [
# 	('query', 'query', 'text'),					# query seqeunce ID
# 	('target', 'target', 'text'),				# matching seqeunce ID
# 	('id', 'identity', 'real'),					# percent identity of the alignment
# 	('ql', 'query_length', 'integer'),			# length of the query sequence
# 	('pairs', 'align_length', 'integer'),		# number of non-gap columns in the alignment
# 	# ('tstrand', 'strand', 'text'),
# 	('tlo', 'match_start', 'integer'), 			# start location in matching sequence
# 	('thi', 'match_end', 'integer'), 			# ending location
# 	('trow', 'target_chars', 'text'),			# sequence chars of match (including gaps)
# ]

# This list has the names of attributes printed to the output file;
# names are passed to the 'userfields' option to usearch

hit_fields = [
	'query',			# query seqeunce ID
	'target', 			# matching seqeunce ID
	'id', 				# percent identity of the alignment
	'ql', 				# length of the query sequence
	'pairs', 			# number of non-gap columns in the alignment
	'tlo',				# start location in matching sequence
	'thi', 				# ending location 
	'trow', 			# sequence chars of match (including gaps)
]

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
# Initialize (or re-initialize) the table that will hold the output from uclust

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)

# Deprecated...
# db_cols = ', '.join(map(lambda x: x[1] + ' ' + x[2], hit_fields))
# create_table = 'CREATE TABLE {tbl} ( hit_id INTEGER PRIMARY KEY AUTOINCREMENT, {cols} )'.format(tbl = table_name, cols = db_cols)

create_table = 'CREATE TABLE {tbl} ( hit_id INTEGER PRIMARY KEY AUTOINCREMENT, panda_id INTEGER, match_id TEXT, identity REAL, query_length INTEGER, match_pairs INTEGER, match_start INTEGER, match_end INTEGER, match_chars TEXT)'.format(tbl = table_name)
	
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
# Print sequences in the uniq table in fasta format.

fetch_sequences = 'SELECT defline, sequence, n FROM uniq JOIN panda USING (panda_id)'
	
def print_sequences(db, args):
	sql = fetch_sequences
	if not args.singletons:
		sql += ' WHERE n > 1'
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', sql))
	
	ff = open(os.path.join(args.workspace, input_file), 'w')
	for defline, sequence, n in db.execute(sql).fetchall():
		print('>{}'.format(defline), file=ff)
		print(sequence, file=ff)
	ff.close()

###
# Run usearch_local, saving results in TSV format.

def run_usearch_local(args):
	cmnd = 'usearch -usearch_local '
	cmnd += os.path.join(args.workspace, input_file)
	cmnd += ' -db ' + args.reference
	cmnd += ' -id {} -strand plus -maxaccepts 0 -maxrejects 64 -maxhits 1'.format(args.identity)
	cmnd += ' -userout ' + os.path.join(args.workspace, output_file)
	# cmnd += ' -userfields ' + '+'.join(map(lambda x: x[0], hit_fields))
	cmnd += ' -userfields ' + '+'.join(hit_fields)
	print(cmnd)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
	db.commit()
	res = os.system(cmnd)
	
###
# Import the TSV file into the hits table.  Easy to use '.import' from the command
# line, can't find how to do it from shell (needs two commands, one to set the 
# separator to tab, second to run the .import command)

# Deprecated...
# cols = ', '.join(map(lambda x: x[1], hit_fields))
# qmarks = ', '.join(map(lambda x: '?',hit_fields))		# baroque way of making one ? for each column....
# insert_hit = 'INSERT INTO {tbl} ({cols}) VALUES ({qmarks})'.format(tbl = table_name, cols = cols, qmarks = qmarks)

insert_hit = 'INSERT INTO {tbl} (panda_id, match_id, identity, query_length, match_pairs, match_start, match_end, match_chars) VALUES (?,?,?,?,?,?,?,?)'.format(tbl = table_name)

def import_results(db, args):
	dm = defline_map(db)					# maps query IDs (deflines) into panda_ids
	qlx = hit_fields.index('ql')			# where to find query length in output record
	pairx = hit_fields.index('pairs')		# where to find number of matching columns
	cutoff = float(args.match)
	for row in open(os.path.join(args.workspace, output_file)):
		cols = row.strip().split('\t')
		if float(cols[pairx]) / float(cols[qlx]) > cutoff:
			db.execute(insert_hit, tuple([dm[cols[0]]] + cols[1:]))
		# else:
		# 	print(cols)

###
# Top level function: initialize the workspace directory, run the app

def form_otus(db, args):
	init_workspace(args)
	print_sequences(db, args)
	run_usearch_local(args)
	import_results(db, args)

###
# Parse the command line arguments, call the top level function...

def init_api():
	parser = argparse.ArgumentParser(
		description="Run usearch to find OTU seeds that match known reference sequences.",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
	parser.add_argument('--singletons', action='store_true', help = 'include singletons (default: disregard singletons)')
	parser.add_argument('-w', '--workspace', help='working directory', default='hits')
	parser.add_argument('-r', '--reference', help='FASTA file containing reference sequences', default=sys.path[0]+'/resources/rRNA16S.gold.fasta')
	parser.add_argument('-i', '--identity', help='minimum percent identity', default=0.97)
	parser.add_argument('-m', '--match', help='minimum percent of query that must align', default=0.95)
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()

	db = sqlite3.connect(args.dbname)
		
	if not prepare_table(db, args):
		argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')

	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	form_otus(db, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

	db.commit()
