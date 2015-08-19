#! /usr/bin/env python3

# Import reads into a SQLite database for a PIP/NGS analysis pipeline.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

# Read a set of fastq files, extract the index from the defline, pack
# the sequence and quality chars into a "BLOB"

import sqlite3
import argparse
import os.path
import sys

from common import *
from FASTQ import *

###
# Make a list of sample IDs to import.  If --sample X was specified on the 
# command line load data for that sample only, otherwise load data for all
# samples.  The return value is a list of records with FASTQ file names for 
# the specified samples.

def sample_list(db, args):
    "Return a list of sample names and the corresponding FASTQ file names currently in the DB"
    sql = 'SELECT sample_id, name, r1_file, r2_file FROM samples'
    if isinstance(args.sample, str):
        sql += ' WHERE name = "{}"'.format(args.sample)
    return db.execute(sql).fetchall()

###
# Prepare 'reads', the table that will hold the FASTQ data.  If the table
# doesn't exist create it.  If we already have a table get the set of sample
# names for reads already in the database and compare them with the names 
# specified by the command line.  If --force was specified on the command
# line erase any conflicing data, otherwise return False so the script prints
# a warning message and exits.

find_read_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "reads"'
create_read_table = 'CREATE TABLE reads ( read_id INTEGER PRIMARY KEY AUTOINCREMENT, sample_id INTEGER NOT NULL REFERENCES samples, read BLOB, defline text DEFAULT NULL )'
find_sample = 'SELECT sample_id FROM reads WHERE sample_id = ? LIMIT 1'
delete_sample = 'DELETE FROM reads WHERE sample_id = ?'
create_read_index = 'CREATE INDEX defx ON reads (defline)'

def prepare_read_table(db, existing, args):
    "Create a new reads table or (if specified by --force) erase previous reads."
    if not db.execute(find_read_table).fetchall():
        db.execute(create_read_table)
        db.execute(create_read_index)
        return True
        
    already_loaded = []
    for x in existing:                # Important: sample ID must be first column in each record
        if db.execute(find_sample, (x[0],)).fetchall():
            already_loaded.append(x[0])
            
    if already_loaded and not args.force:
        return False
        
    for sample_id in already_loaded:        # if none previousy loaded this is a NOP
        db.execute(delete_sample, (sample_id,)) 
           
    return True     

###
# Read the sequences from a pair of files (R1 and R2) for a specified sample.  Create
# a FASTQ object for each read, compress it, and insert the compressed form into the 
# reads table.

insert_blob = 'INSERT INTO reads (sample_id, read) VALUES (?, ?)'
insert_blob_and_def = 'INSERT INTO reads (sample_id, read, defline) VALUES (?, ?, ?)'
    
def load_sequences(db, fn1, fn2, sid, args):
    "Load files fn1 and fn2 into the reads table"
    
    if not (os.path.exists(fn1) and os.path.exists(fn2)):
        return
    
    file1 = FASTQReader(fn1)
    file2 = FASTQReader(fn2)
    count = 0
    
    limit = int(args.limit) if args.limit else None

    seq1 = file1.readseq()
    seq2 = file2.readseq()
    
    while seq1 is not None and seq2 is not None:
        for seq in [seq1, seq2]:
            if args.quality and seq.filtered() == 'Y':
                blob = None
            else:
                seq.pack()
                blob = seq.blob()
            if args.deflines:
                # db.execute(insert_blob_and_def, (sid, blob, seq.unique_defline()))
                db.execute(insert_blob_and_def, (sid, blob, seq.defline()))
            else:
                db.execute(insert_blob, (sid, blob))
        count += 1
        if limit and count >= limit:
            break
        seq1 = file1.readseq()
        seq2 = file2.readseq()
        
    file1.close()
    file2.close()
    
###
# The main function -- iterate over the names of the samples, call
# load_sequences to import the sequences from the FASTQ files.
    
def import_files(db, samples, args):
    for sid, sname, fastq1, fastq2 in samples:
        fn1 = os.path.join(args.directory, fastq1)
        fn2 = os.path.join(args.directory, fastq2)
        if not args.noimport:
            load_sequences(db, fn1, fn2, sid, args)
            record_metadata(db, 'import', '{}, {}'.format(fn1,fn2))

###
# Check the combination of command line options to make sure they're sensible

def validate_options(args):
    # if we're not loading data the limit, defline, and quality options are superfluous
    if args.noimport and (args.limit or args.deflines or args.quality):
        print('Warning: options ignored: with --noimport the following options are ignored: --limit, --defline, --quality')

###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Import FASTQ files into a SQLite3 database for a PIP/NGS analysis pipeline.",
        epi = "If --noimport is specified sequence descriptions are created but sequences are not loaded into the database.",
        specs = [
            ('directory',    { 'required' : True, 'metavar': 'dir', 'help' : '(required) name of directory containing FASTQ files' } ),
            ('noimport', { 'action': 'store_true', 'help' : "record file names but don't import sequences"} ),
            ('sample', { 'metavar': 'id', 'help' : 'import sequences for this sample only'} ),
            ('update', { 'action': 'store_true', 'help' : 'used with --sample, add sequences to existing table'} ),
            ('limit', { 'metavar': 'N', 'help' : 'max number of sequences imported for each sample'} ),
            ('quality', { 'action': 'store_true', 'help' : "don't import sequences flagged as low quality"} ),
            ('deflines', { 'action': 'store_true', 'help' : 'include deflines with sequences'} ),
        ]
    )
    validate_options(args)
    
    db = sqlite3.connect(args.dbname)
    
    samples = sample_list(db, args)
    if not prepare_read_table(db, samples, args):
        argparse.ArgumentParser.exit(1, 'Sequences were loaded previously; use --force to replace them')
                
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))
    import_files(db, samples, args)
    record_metadata(db, 'end', '')
    
    db.commit()
