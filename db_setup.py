#! /usr/bin/env python3

# Load barcodes and other auxiliary info into a SQLite database for a PIP/NGS analysis pipeline.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

# Usage:
#
#    db_setup.py [--force] dbname --info experiment
#
# Add sequence information to the database file with the specified name.  The --info argument
# is the name of a file with barcodes, index sequences, and other information 
# about the sequences.  By default the script checks to make sure the database 
# does not exist already; use --force to re-initialize an existing database.

import sqlite3
import argparse
import os.path
import sys

# Strings in this list are names of items from the info file; plural forms
# will be names of tables in the database

items = ['primer', 'offset']

# SQL command templates

create_table = "CREATE TABLE {}s ( pair_end INTEGER, sequence TEXT )"
drop_table = "DROP TABLE IF EXISTS {}s"
insert = "INSERT INTO {}s VALUES (?, ?)"

### 
# Helpers

def init_log(db):
    "Create the log table, insert initial log messages"
    db.execute("DROP TABLE IF EXISTS log")
    db.execute("CREATE TABLE log ( time timestamp, script text, event text, message text )")
    
def make_table(db, tbl):
    "Initialize a table -- delete it if it exists already"
    db.execute(drop_table.format(tbl))
    db.execute(create_table.format(tbl))
    
def add_aux_info(db, args):
    "Create primer and offset tables, read the info file and populate the tables"
    
    for tbl in items:
        make_table(db, tbl)
    
    if args.info is not None:
        for line in open(args.info):
            tbl, loc, seq = line.strip().split()
            if tbl in ['primer', 'offset']:
                db.execute('INSERT INTO {}s VALUES (?, ?)'.format(tbl), (loc, seq))
            else:
                print("Unknown item:", line.strip())
            
def add_sample_info(db, args):
    "Create the sample table with sample names and the names of their FASTQ files"
    
    db.execute('DROP TABLE IF EXISTS samples')
    db.execute('CREATE TABLE samples (sample_id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, r1_file TEXT, r2_file TEXT, barcode TEXT )')

    for line in open(args.data):
        a = line.strip().split()
        if len(a) < 4:
            a.append(None)
        db.execute('INSERT INTO samples (name, r1_file, r2_file, barcode) VALUES (?, ?, ?, ?)', tuple(a))

### 
# Set up command line arguments

def init_api():
    parser = argparse.ArgumentParser(
        description="Initialize a SQLite3 database for a PIP/NGS analysis pipeline.",
        epilog="If --force is not specified the program will print a message and exit if the database exists already."
    )
    parser.add_argument('dbname', help='the name of the SQLite database file')
    parser.add_argument('--info', help='[GUI:path:info:<req>] path to file with primer and offset specs')
    parser.add_argument('--data', required=True, help='[GUI:path:data:<req>] path to file with sample names and corresponding FASTQ file names')
    parser.add_argument('--force', action='store_true', help='[GUI:flag::T] re-initialize an existing database')
    parser.add_argument('--message', metavar='text', help='[GUI:str:Project description:] initial log message')
    return parser.parse_args()
    
###
# Top level....
    
if __name__ == "__main__":
    args = init_api()

    if os.path.exists(args.dbname) and not args.force:
        argparse.ArgumentParser.exit(1, 'Found existing database with that name; use --force to reinitialize')

    db = sqlite3.connect(args.dbname)
    
    init_log(db)
    add_aux_info(db, args)
    add_sample_info(db, args)
    
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
    if args.message is not None:
        db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'message', args.message))

    db.commit()
