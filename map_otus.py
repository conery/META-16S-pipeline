#! /usr/bin/env python3

# Run usearch to map sequences to OTUs

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

from FASTQ import *

import sqlite3
import argparse
import os
import os.path
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'otu'
result_file_pattern = 'readmap.{}.uc'
input_file_pattern = 'unique.{}.fasta'
ref_db_file = 'otus.fasta'

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
create_table = 'CREATE TABLE {tbl} ( otu_id INTEGER, sample_id INTEGER, count INTEGER )'.format(tbl = table_name)
	
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
# Make the reference "database" (FASTQ file) from non-chimeric OTU seeds

# fetch_sequences = 'SELECT otu_id, sequence FROM seeds JOIN chimeras USING (panda_id) WHERE chimeric = "N"'
fetch_sequences = "SELECT cluster_id, sequence FROM (SELECT cluster_id, count(members.name) AS n, sequence FROM clusters LEFT JOIN members USING (cluster_id) GROUP BY cluster_id) JOIN chimeras USING (cluster_id) WHERE n > 0 AND chimeric = 'N'"

def make_reference_db(db, args):
	ff = open(os.path.join(args.workspace, ref_db_file), 'w')
	for otu_id, sequence in db.execute(fetch_sequences).fetchall():
		print('>OTU_{}'.format(otu_id), file=ff)
		print(sequence, file=ff)
	ff.close()

###
# Run the app

def run_usearch_global(sid, args):
	cmnd = 'usearch -usearch_global '
	cmnd += os.path.join(args.directory, input_file_pattern.format(sid))
	cmnd += ' -db ' + os.path.join(args.workspace, ref_db_file)
	cmnd += ' -strand plus'
	cmnd += ' -id 0.97'
	cmnd += ' -uc ' + os.path.join(args.workspace, result_file_pattern.format(sid))
	print(cmnd)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
	db.commit()
	res = os.system(cmnd)

###
# Fetch a list of sample ids 

# fetch_barcodes = 'SELECT barcode_id, experiment FROM barcodes'

# def load_barcodes(db, args):
# 	files = [ ]
# 	sql = fetch_barcodes
# 	if isinstance(args.experiment, str):
# 		sql += ' WHERE experiment = "{}"'.format(args.experiment)
# 	return list(map(lambda x: x[0], db.execute(sql).fetchall()))
	
fetch_samples = 'SELECT sample_id FROM samples'

def fetch_sample_list(db, args):
	"Return a list of sample IDs for the samples to process"
	sql = fetch_samples
	if isinstance(args.sample, str):
		sql += ' WHERE name = "{}"'.format(args.sample)
	return db.execute(sql).fetchall()


# Populate the table 

insert_record = 'INSERT INTO {} (otu_id, sample_id, count) VALUES (?,?,?)'.format(table_name)
fetch_counts = 'SELECT defline, n FROM panda JOIN uniq USING (panda_id) where uniq.sample_id = {}'

def import_results(db, args, sid):
	cmap = defline_map(db, sid)			# map deflines to number of times seq found in this sample
	count = { }
	for line in open(os.path.join(args.workspace, result_file_pattern.format(sid))):
		res = line.split('\t')
		otu = res[-1].strip()
		otu_id = 0 if otu == '*' else int(otu.split('_')[-1])
		defline = FASTQ.parse_defline(res[8])
		count.setdefault(otu_id, 0)
		count[otu_id] += cmap[defline]
	for otu_id in sorted(count.keys()):
		db.execute(insert_record, (otu_id, sid, count[otu_id]))

def defline_map(db, sid):
	dm = { }
	query = fetch_counts.format(sid)
	for defline, n in db.execute(query):
		dm[defline] = n
	return dm
	
###
# Top level function: initialize the workspace directory, run the app

def map_otus(db, args):
	init_workspace(args)
	make_reference_db(db,args)
	for row in fetch_sample_list(db, args):
		sid = row[0]
		run_usearch_global(sid, args)
		import_results(db, args, sid)

###
# Parse the command line arguments, call the top level function...

def init_api():
	parser = argparse.ArgumentParser(
		description="Run usearch to map merged sequences to one of the inferred clusters",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
	parser.add_argument('-w', '--workspace', help='working directory', default='map')
	parser.add_argument('-d', '--directory', help='directory with input file', default='uniq')
	parser.add_argument('-s', '--sample', metavar='id', required=False, help='map only this sample')
	# parser.add_argument('-r', '--reference', help='FASTA file containing reference sequences', default='reference/gold.fa')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()

	db = sqlite3.connect(args.dbname)
		
	if not prepare_table(db, args):
		argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')

	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	map_otus(db, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

	db.commit()
