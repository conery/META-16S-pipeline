#! /usr/bin/env python3

# Import reads into a SQLite database for a PIP/NGS analysis pipeline.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

# Read a set of fastq files, extract the index from the defline, pack
# the sequence and quality chars into a "BLOB"

from common import *

import sqlite3
import argparse
import os.path
import sys

from FASTQ import *

fetch_samples = 'SELECT sample_id, name, r1_file, r2_file FROM samples'

def sample_list(db, args):
	"Return a list of sample names and the corresponding FASTQ file names."
	sql = fetch_samples
	if isinstance(args.sample, str):
		sql += ' WHERE name = "{}"'.format(args.sample)
	return db.execute(sql).fetchall()
	
find_read_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "reads"'
create_read_table = 'CREATE TABLE reads ( read_id INTEGER PRIMARY KEY AUTOINCREMENT, sample_id INTEGER NOT NULL REFERENCES samples, read BLOB, defline text DEFAULT NULL )'
find_sample = 'SELECT sample_id FROM reads WHERE sample_id = ? LIMIT 1'
delete_sample = 'DELETE FROM reads WHERE sample_id = ?'
create_read_index = 'CREATE INDEX defx ON reads (defline)'

def prepare_tables(db, samples, args):
	"""
	Return True if the database is initialized and ready to accept a new set of reads.
	If the reads table does not exist yet create it and return True.  If the 
	samples have already been loaded return False unless --force was specified on the
	command line.  If --force was specified erase the old data and return True.
	"""
	if not db.execute(find_read_table).fetchall():
		db.execute(create_read_table)
		db.execute(create_read_index)
		return True
	already_loaded = []
	for x in samples:				# Important: sample ID must be first column in each record
		if db.execute(find_sample, (x[0],)).fetchall():
			already_loaded.append(x[0])
	if already_loaded and not args.force:
		return False
	for sid in already_loaded:		# if none previousy loaded this is a NOP
		db.execute(delete_sample, (sid,))	
	return True	 
	
insert_blob = 'INSERT INTO reads (sample_id, read) VALUES (?, ?)'
insert_blob_and_def = 'INSERT INTO reads (sample_id, read, defline) VALUES (?, ?, ?)'
	
def load_sequences(db, fn1, fn2, sid, args):
	"""
	Read the sequences from a pair of files (R1 and R2) for a specified sample.  Create
	a FASTQ object for each read, compress it, and insert the compressed form into the 
	reads table.
	
	TBD: if the --accelerate option is specified on the command line use the C++ app
	to fill the table.
	"""
	if not (os.path.exists(fn1) and os.path.exists(fn2)):
		return
	
	# print(fn1)
	
	file1 = open(fn1)
	file2 = open(fn2)
	count = 0
	
	limit = int(args.limit) if args.limit else None

	seq1 = FASTQ.read(file1)
	seq2 = FASTQ.read(file2)
	while seq1 is not None and seq2 is not None:
		for seq in [seq1, seq2]:
			# defline = seq.defline()
			# if args.quality and defline.find('Y', defline.index(' ')) > 0:
			if args.quality and seq.filtered() == 'Y':
				blob = None
			else:
				seq.pack()
				blob = seq.blob()
			if args.defline:
				# seq_id, flags = defline.split()
				db.execute(insert_blob_and_def, (sid, blob, seq.unique_defline()))
			else:
				db.execute(insert_blob, (sid, blob))
		count += 1
		if limit and count >= limit:
			break
		seq1 = FASTQ.read(file1)
		seq2 = FASTQ.read(file2)
		
	file1.close()
	file2.close()
	
# file_pattern = 'L1_{exp}_{code1}-{code2}_L001_R{pair}_001.fastq'
	
def import_files(db, samples, args):
	"""
	Use the barcodes to create the names of the files to load, pass the file names to
	load_sequences to read the sequences and save them in the database.
	"""
	for sid, sname, fastq1, fastq2 in samples:
		fn1 = os.path.join(args.directory, fastq1)
		fn2 = os.path.join(args.directory, fastq2)
		if not args.noimport:
			load_sequences(db, fn1, fn2, sid, args)
		db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'paired', '{}, {}'.format(fn1,fn2)))
		
### 
# Set up command line arguments

def init_api():
	parser = argparse.ArgumentParser(
		description="""Import FASTQ files into a SQLite3 database for a PIP/NGS analysis pipeline.  If the
		--noimport option is specified only descriptions of the FASTQ files are saved.
		""",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('--directory', metavar='dir', required=True, help='[GUI:dir:Directory containing FASTQ files:<req>] name of data directory')
	parser.add_argument('--limit', metavar='N', required=False, type=int, help='[GUI:int:Limit:] maximum number of reads per experiment')
	parser.add_argument('--sample', metavar='id', required=False, help='import files for sample id')
	parser.add_argument('--defline', action='store_true', help='[GUI:flag:Save deflines?:F] include deflines')
	parser.add_argument('--quality', action='store_true', help="[GUI:flag:High quality only?:F] filter out low quality reads")
	parser.add_argument('--force', action='store_true', help='[GUI:flag::T] re-initialize an existing database')
	parser.add_argument('--all', action='store_true', help='[GUI:flag::T] import files for all experiments')
	parser.add_argument('--noimport', action='store_true', help="record file names but don't import")
	parser.add_argument('--warn', action='store_true', help="print a warning message for files not loaded")
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	
	if not (args.all or args.sample):
		argparse.ArgumentParser.exit(1, 'Specify a sample name with --sample ID or import all with --all')

	db = sqlite3.connect(args.dbname)
	
	samples = sample_list(db, args)
	
	if not prepare_tables(db, samples, args):
		argparse.ArgumentParser.exit(1, 'These sequences were loaded previously; use --force if you want to replace them')
				
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	import_files(db, samples, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )
	
	db.commit()
