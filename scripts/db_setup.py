#! /usr/bin/env python3

# Create a new SQLite database for a PIP/NGS analysis pipeline.

# John Conery / Kevin Xu Junjie
# University of Oregon

# Usage:
#
#    db_setup.py dbname [--force] --data samples.csv [--primers primers.csv]
#
# The required arguments are the database name and --data, which specifies a CSV
# file that contains the names of the samples and the FASTQ files containing the
# data for each sample.  
#
# If the data includes primers or offsets they can be specified in a second TSV file.  
#
# Use --force to re-initialize an existing database.
#
# Another optional argument, --message, can be used to add a project description
# message to the log table.

import sqlite3
import argparse
import os.path
import sys

from common import *

### 
# Create the metadata table

def init_log(db):
    "Create the log table, insert initial log messages"
    db.execute("DROP TABLE IF EXISTS log")
    db.execute("CREATE TABLE log ( time timestamp, script text, event text, message text )")
    
###
# Initialize the primers and offsets tables.  If the user specified a TSV file
# with data for these tables add the records, otherwise leave the tables empty. 
    
create_table = "CREATE TABLE {}s ( pair_end INTEGER, sequence TEXT )"
drop_table = "DROP TABLE IF EXISTS {}s"

def make_table(db, tbl):
    "Initialize a table -- delete it if it exists already"
    db.execute(drop_table.format(tbl))
    db.execute(create_table.format(tbl))
    
def add_aux_info(db, args):
    "Create primer and offset tables, read the info file and populate the tables"
    
    for tbl in ['primer', 'offset']:
        make_table(db, tbl)
    
    if args.primers is not None:
        for line in open(args.primers):
            tbl, loc, seq = line.strip().split()
            if tbl in ['primer', 'offset']:
                db.execute('INSERT INTO {}s VALUES (?, ?)'.format(tbl), (loc, seq))
            else:
                print("Unknown item:", line.strip())

###
# Create the table that describes samples in the experiment.  Each record has
# a sample ID and the names of the two FASTQ files with data for that sample.
            
def add_sample_info(db, args):
    "Create the sample table with sample names and the names of their FASTQ files"
    
    db.execute('DROP TABLE IF EXISTS samples')
    db.execute('CREATE TABLE samples (sample_id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, r1_file TEXT, r2_file TEXT)')

    for line in open(args.data):
        a = line.strip().split()
        db.execute('INSERT INTO samples (name, r1_file, r2_file) VALUES (?, ?, ?)', tuple(a))


###
# Top level....
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Initialize a SQLite3 database for a PIP/NGS analysis pipeline.",
        epi = None,
        specs = [
            ('data',    { 'required' : True, 'metavar': 'x.csv', 'help' : '(required) CSV file with names of samples and their FASTQ files' } ),
            ('primers', { 'metavar': 'x.csv', 'help' : '(optional) CSV file with primer and/or offset sequences' } ),
            ('message', { 'metavar': '"text"', 'help' : 'project description to add to metadata'} ),
        ]
    )

    if os.path.exists(args.dbname) and not args.force:
        argparse.ArgumentParser.exit(1, 'Found existing database with that name; use --force to reinitialize')
        
    db = sqlite3.connect(args.dbname)
    
    init_log(db)
    add_aux_info(db, args)
    add_sample_info(db, args)
    
    record_metadata(db, 'init', ' '.join(['data:', args.data, 'primers:', args.primers or 'None']))
    record_metadata(db, 'desc', args.message)

    db.commit()
