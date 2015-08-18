#! /usr/bin/env python3

# Run cd-hit-dup to find unique assembled reads.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-05-07

from common import *
from FASTQ import *

import sqlite3
import argparse
import os
import os.path
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'uniq'
temp_table_name = 'clstr'
seq_table_name = 'panda'

input_file_pattern = 'merged.{}.fasta'
output_file_pattern = 'unique.{}.fasta'

###
# Utility -- echo a query on the terminal, add it to the log table, and run it

def print_log_and_query(cmnd):
    print(cmnd)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', cmnd))
    db.execute(cmnd)
    
###
# Initialize (or re-initialize) the table that will hold the output from cd-hit-dup
    
find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)
clear_seq_table = 'DELETE FROM {}'.format(seq_table_name)
# create_temp_table = 'CREATE TEMPORARY TABLE {tbl} ( {tbl}_id INTEGER PRIMARY KEY AUTOINCREMENT, bcid INTEGER, defline TEXT, n INTEGER )'.format(tbl = temp_table_name)
create_table = 'CREATE TABLE {tbl} ( {tbl}_id INTEGER PRIMARY KEY AUTOINCREMENT, panda_id INTEGER NOT NULL REFERENCES panda, sample_id INTEGER, n INTEGER )'.format(tbl = table_name)

def prepare_tables(db, args):
    """
    Check to see if the DB already has a 'uniq' table, and return False unless the user
    put --force on the command line.  Also, clear the sequence table so we can import
    new ones if the --load_seqs flag is set.
    """
    if db.execute(find_table).fetchall():
        if args.force:
            db.execute(drop_table)
        else:
            return False
    db.execute(create_table)
    if args.load_seqs:
        db.execute(clear_seq_table)
    return True
    
    
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
# Fetch a list of sample ids (will be used to format file names)

# fetch_barcodes = 'SELECT barcode_id, experiment FROM barcodes'
#
# def load_barcodes(db, args):
#     files = [ ]
#     sql = fetch_barcodes
#     if isinstance(args.experiment, str):
#         sql += ' WHERE experiment = "{}"'.format(args.experiment)
#     return list(map(lambda x: x[0], db.execute(sql).fetchall()))

fetch_samples = 'SELECT sample_id FROM samples'

def fetch_sample_list(db, args):
    "Return a list of sample IDs for the samples to process"
    sql = fetch_samples
    if isinstance(args.sample, str):
        sql += ' WHERE name = "{}"'.format(args.sample)
    return db.execute(sql).fetchall()

###
# Run cd-hit-dup using options specified on the command line
# TBD: write to /dev/null?  -N?? 

def run_cd_hit_dup(sid, args):
    cmnd = 'cd-hit-dup'
    cmnd += ' -i ' + os.path.join(args.directory, input_file_pattern.format(sid))
    cmnd += ' -o ' + os.path.join(args.workspace, output_file_pattern.format(sid))
    print(cmnd)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
    db.commit()
    res = os.system(cmnd)
    
###
# Open the cluster file, get info about each cluster, return as a dictionary.
# The helper 'get_cluster' reads the lines in one cluster, returning the 
# line that follows a cluster and the lines in a cluster.  Save cluster info
# in the dictionary.

def get_cluster(line, file, info):
    cid = int(line.split()[1])        # first line is '>Cluster N'
    line = file.readline()
    seqinfo = line.split()[2]        # next line has defline of first seq in the cluster
    # defline = ':'.join(seqinfo[1:].split(':')[:-1])
    defline = FASTQ.parse_defline(seqinfo)
    count = 0
    while len(line) > 0 and line[0] != '>':
        count += 1
        line = file.readline()
    info[defline] = count
    return line

def parse_clusters(args, sid):
    clusters = { }
    file = open(os.path.join(args.workspace, output_file_pattern.format(sid) + '.clstr'))
    line = file.readline()
    while len(line) > 0:
        line = get_cluster(line, file, clusters)
    return clusters


###
# Load sequences from the fasta file produced by cd-hit-dup into the panda (paired
# sequence) table

insert_sequence = 'INSERT INTO {} (sample_id, defline, sequence) VALUES (?, ?, ?)'.format(seq_table_name)

def import_sequences(db, args, sid):
    file = open(os.path.join(args.workspace, output_file_pattern.format(sid)))
    defline = file.readline()
    while len(defline) > 0:
        seq = file.readline()
        # strip barcode sequence from defline
        # parts = defline.split(':')
        # defline = ':'.join(parts[:7])
        defline = FASTQ.parse_defline(defline)
        db.execute(insert_sequence, (sid, defline, seq.strip()))
        defline = file.readline()

###
# Two things to do after each run of cd-hit-dup: if the panda table is empty we
# need to import the sequences, then import the clusters into the uniq table.

# insert_cluster = 'INSERT INTO {} (bcid, defline, n) VALUES (?, ?, ?)'.format(temp_table_name)
insert_cluster = 'INSERT INTO {} (panda_id, sample_id, n) VALUES (?, ?, ?)'.format(table_name)

def import_results(db, args, sid):
    if args.load_seqs:
        import_sequences(db, args, sid)
    defmap = defline_map(db, sid)
    clusters = parse_clusters(args, sid)
    count = 0
    limit = int(args.limit) if args.limit else None
    for defline, size in clusters.items():
        # db.execute(insert_cluster, (bcid, defline, size))
        db.execute(insert_cluster, (defmap[defline], sid, size))
        count += 1
        if limit and count >= limit:
            break

###
#     --- old description ---
# Top level function: initialize the workspace directory, figure out which barcodes
# to use and the number of sequences already in the panda table, and then run the
# app for each barcode.  After each run cluster data is loaded into a temp table.
# The final step creates the output table by combining columns from the temp table
# with sequence data.

# fix -- create-select, use global var table names
# create_table = 'CREATE TABLE {tbl} ( {tbl}_id INTEGER PRIMARY KEY AUTOINCREMENT, panda_id INTEGER NOT NULL REFERENCES panda, n INTEGER )'.format(tbl = table_name)

# create_index = 'CREATE INDEX {idx} on {tbl} (defline)'
# create_table = 'CREATE TABLE {tbl} AS SELECT {tmp_tbl}_id AS {tbl}_id, {seq_tbl}_id, barcode_id, n FROM {tmp_tbl} JOIN {seq_tbl} USING (defline)'.format(tbl = table_name, tmp_tbl = temp_table_name, seq_tbl = seq_table_name)

###
# Top level function: run cd-hit-dup for each set of paired sequences, then save
# the id of the first sequence in the cluster and the number of sequences in the
# cluster.

def remove_duplicates(db, args):
    init_workspace(args)
    for row in fetch_sample_list(db, args):
        sid = row[0]
        run_cd_hit_dup(sid, args)
        import_results(db, args, sid)
        
    # print_log_and_query(create_index.format(idx = 'pdef', tbl = 'panda'))
    # print_log_and_query(create_index.format(idx = 'cdef', tbl = 'clstr'))
    # print_log_and_query(create_table)

### 
# Set up command line arguments

def init_api():
    parser = argparse.ArgumentParser(
        description="""Run cd-hit-dup to find unique assembled sequences.
        """
    )
    parser.add_argument('dbname', help='the name of the SQLite database file')
    parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
    parser.add_argument('--load_seqs', action='store_true', help='populate the panda table with uniq sequences')
    parser.add_argument('-w', '--workspace', help='working directory', default='uniq')
    parser.add_argument('-d', '--directory', help='directory with input FASTQ files', default='panda')
    parser.add_argument('-s', '--sample', metavar='id', required=False, help='process sequences from this sample only')
    parser.add_argument('-l', '--limit', help='number of sequences to load')
    return parser.parse_args()
    
###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    args = init_api()

    db = sqlite3.connect(args.dbname)
    
    q = 'SELECT count(*) FROM {}'.format(seq_table_name)
    seqs_loaded = db.execute(q).fetchall()[0][0]
    if seqs_loaded == 0 and not args.load_seqs:
        argparse.ArgumentParser.exit(1, 'No paired sequences in panda table; specify --load_seqs or re-run assemble_pairs without --noimport')
        
    if not prepare_tables(db, args):
        argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')
        
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
    remove_duplicates(db, args)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )
        
    db.commit()
